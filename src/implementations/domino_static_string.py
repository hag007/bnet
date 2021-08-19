import sys

sys.path.insert(0, "../../")

import random
import os
import pandas as pd
import numpy as np
import pickle
import multiprocessing
from functools import reduce
from networkx.algorithms.boundary import edge_boundary
from scipy.stats import hypergeom
from statsmodels.sandbox.stats.multicomp import fdrcorrection0

import pcst_fast
import networkx as nx
from networkx.algorithms.community.quality import modularity
from networkx.algorithms.community.centrality import girvan_newman
from networkx.algorithms.components import connected_components

from src.utils.graph_influence_linear_th import linear_threshold
from src.implementations.preprocess_slices import read_preprocessed_slices
from src.implementations.network_builder import build_network_static

import src.constants as constants

G_modularity = None

ACTIVITY_BASELINE=0 # 0.2
SIMILARITY_FACTOR=2.5 # 0.2
def extract_scores(scores_file):
    """"""
    scores = pd.read_csv(scores_file, sep='\t', index_col=0)
    if "qval" in scores.columns:

        scores["score"] = -np.log10(scores["qval"])/1.3
    else:
        scores = pd.read_csv(scores_file, sep='\t', index_col=0, header=None, dtype=str)
        scores["score"] = 1
    return scores


def add_scores_to_nodes(G, scores):
    """"""
    inds = []
    for nd in G.nodes:
        G.nodes[nd]["pertubed_node"] = False
        G.nodes[nd]["score"] = 0

    c=0
    for ind, row in scores.iterrows():
        if ind in G.nodes:
            inds.append(ind)
            G.nodes[ind]["score"] = row["score"]
            G.nodes[ind]["pertubed_node"] = (row["score"] >= 1)  # binarizing the activeness
            c+=int(row["score"] >= 1)

    print(f'c={c}')

    return G


def create_subgraph(params):
    cur_module = params
    global G_modularity
    nodes = set(cur_module)
    res = G_modularity.subgraph(list(nodes))
    return res


def prune_network_by_modularity(G, modules, cache_file):
    global G_modularity
    if os.path.exists(cache_file) and constants.USE_CACHE:
        print(f'fetch cache file for subnetworks {cache_file}')
        G_modularity = pickle.load(open(cache_file, 'rb'))
        for n in G_modularity:
            G_modularity.nodes[n]['pertubed_node'] = G.nodes[n]['pertubed_node']
        print('pkl is loaded')
        return G_modularity

    print(f'generating subgraphs...')
    G_modularity = G
    print(f"Before slicing: n of cc:{len(list(connected_components(G_modularity)))}, n of nodes: {len(G_modularity.nodes)}, n of edges, {len(G_modularity.edges)}")
    p = multiprocessing.Pool(constants.N_OF_THREADS)

    G_modules = p.map(create_subgraph, [m for m in modules])
    p.close()
    # print(f'{modules}')
    print(f'# of modules after extraction: {len(G_modules)}')
    G_modularity = nx.algorithms.operators.union_all(G_modules)
    print(f"After slicing: n of cc:{len(list(connected_components(G_modularity)))}, n of nodes: {len(G_modularity.nodes)}, n of edges, {len(G_modularity.edges)}")
    pickle.dump(G_modularity, open(cache_file, 'wb+'))
    print('subgraphs\' pkl is saved')


def prune_network_by_modularity_old(G, modules, dummy):
    G_modularity = G.copy()
    edges_to_remove = []
    for cur_edge in G_modularity.copy().edges:
        in_cc = False
        for cur_module in modules:
            if cur_edge[0] in cur_module and cur_edge[1] in cur_module:
                in_cc = True
        if not in_cc:
            edges_to_remove.append(cur_edge)

    G_modularity.remove_edges_from(edges_to_remove)
    return G_modularity


def get_seeds(slice, threhold):
    get_seeds()


def sa_for_slice(params):
    G, slice, i_slice, ts, relevant_slices, prize_factor, module_threshold = params
    score_sum=sum([G.nodes[a]["score"] for a in slice.nodes])
    subslices = [(G.subgraph(a), G.nodes[a]["score"]/score_sum) for a in slice.nodes]
    for t in np.arange(ts+10, 10,-1):
        random.shuffle(subslices)
        rand_num = np.random.uniform(0, 1, 1)
        prob_agg=0
        for subslice_index, subslice in enumerate(subslices):
            prob_agg+=subslice[1]
            # print(f'cur prob: {prob_agg}')
            if prob_agg >= rand_num:
                updated_subslice_subgraph = modify_subslice(subslices[subslice_index][0], slice, t/40.0)
                break
        subslices=update_subslice_probabilities([a[0] for a in subslices], updated_subslice_subgraph, subslice_index, score_sum)
        # print("current iteration:")
        # print((sum([len(a[0].nodes) for a in subslices]),[len(a[0].nodes) for a in subslices]))

    if sum([len([b for a in subslices for b in a[0].nodes  if slice.nodes[b]["pertubed_node"]])]) > len([a for a in slice if slice.nodes[a]["pertubed_node"]]):
        x=1

    print('subslices fraction for original slice: {}'.format((sum([len(a[0].nodes) for a in subslices]),
    len(slice.nodes),
    sum([len([b for a in subslices for b in a[0].nodes  if slice.nodes[b]["pertubed_node"]])]),
    len([a for a in slice if slice.nodes[a]["pertubed_node"]]),
    [len(a[0].nodes) for a in subslices])))
    return [a[0] for a in subslices]


def update_subslice_probabilities(subslices, new_subslice, modified_subslice_index, score_sum):
    new_subslice_subgraph_list=[]

    for subslice_index, a in enumerate(subslices):
        if modified_subslice_index==subslice_index:
            if len(subslices[subslice_index].nodes) > len(new_subslice):
                removed_node=list(set(subslices[subslice_index].nodes).difference(set(new_subslice)))[0]
                removed_node_G=nx.Graph()
                removed_node_G.add_node(removed_node)
                removed_node_G.nodes[removed_node]["score"]=subslices[subslice_index].nodes[removed_node]["score"]
                removed_node_G.nodes[removed_node]["pertubed_node"] =subslices[subslice_index].nodes[removed_node]["pertubed_node"]
                new_subslice_subgraph_list.append(removed_node_G)

        elif len(set(a.nodes).intersection(new_subslice.nodes))>0 :
            new_subslice=nx.compose(new_subslice,a)
        else:
            new_subslice_subgraph_list.append(a)

    new_subslice_subgraph_list.append(new_subslice)
    return [(cur_subgraph, sum([cur_subgraph.nodes[a]["score"] for a in cur_subgraph.nodes])/score_sum) for cur_subgraph in new_subslice_subgraph_list]



def modify_subslice(cur_subslice, cur_slice, t):
    step = np.random.binomial(1, 0.5)
    modified_slice = None
    if (step or len(cur_subslice.nodes) < 4) and len(cur_subslice.nodes) < len(cur_slice.nodes) :
        modified_slice = add_to_subslice(cur_subslice, cur_slice, t)
    else:
        modified_slice = remove_from_subslice(cur_subslice, t)

    return modified_slice

def get_edge_weight(edge):

    open("/specific/netapp5/gaga/hagailevi/omics/output/weights/static_weights.txt", 'a+').write(
        f'{1/(1+G_modularity.get_edge_data(edge[0],edge[1])["score"]/200.0)}\n')
    return  1/(1+G_modularity.get_edge_data(edge[0],edge[1])["score"]/200.0)

def add_to_subslice(cur_subslice, cur_slice, t):
    cur_slice_minus_cur_subslice=cur_slice.copy()
    cur_slice_minus_cur_subslice.remove_nodes_from(cur_subslice.nodes)
    cross_edges = list(edge_boundary(cur_slice, cur_subslice, cur_slice_minus_cur_subslice))
    edge_to_add = random.choice(cross_edges)
    node_outside_sublice = edge_to_add[0] if edge_to_add[0] not in cur_subslice.nodes else edge_to_add[1]
    # if len([a for a in cur_subslice.nodes if cur_slice.nodes[a]["pertubed_node"]]) ==0:
    #     x=1
    delta_c = -np.mean(get_edge_weight(edge_to_add)) * SIMILARITY_FACTOR + ACTIVITY_BASELINE + \
              float(cur_slice.nodes[node_outside_sublice]['score'])
    # print(f"delta_c addition: {delta_c}")
    p = None
    if delta_c > 0:
        p = 1
    else:
        p = np.e ** (delta_c / t)
    step = np.random.binomial(1,p)
    if step:
        new_subgraph=cur_subslice.copy()
        new_subgraph.add_edge(edge_to_add[0], edge_to_add[1])
        new_subgraph.nodes[edge_to_add[0]]["score"]=cur_slice.nodes[edge_to_add[0]]["score"]
        new_subgraph.nodes[edge_to_add[1]]["score"]=cur_slice.nodes[edge_to_add[1]]["score"]
        new_subgraph.nodes[edge_to_add[0]]["pertubed_node"] = cur_slice.nodes[edge_to_add[0]]["pertubed_node"]
        new_subgraph.nodes[edge_to_add[1]]["pertubed_node"] = cur_slice.nodes[edge_to_add[1]]["pertubed_node"]
        return new_subgraph
    else:
        return cur_subslice.copy()


def remove_from_subslice(cur_subslice, t):
    n_pertubed_nodes=len([a for a in cur_subslice if cur_subslice.nodes[a]["pertubed_node"]])
    leafs = [a for a in cur_subslice if cur_subslice.degree(a) <= 1 and (n_pertubed_nodes>1 or not cur_subslice.nodes[a]["pertubed_node"])]
    if len(leafs)==0:
        return cur_subslice.copy()
    leaf_to_remove = random.choice(leafs)
    # print(cur_subslice.edges(leaf_to_remove))
    edge_to_remove=list(cur_subslice.edges(leaf_to_remove))[0]
    # print(edge_to_remove)
    delta_c = np.mean(get_edge_weight(edge_to_remove)) * SIMILARITY_FACTOR - ACTIVITY_BASELINE - \
              cur_subslice.nodes[leaf_to_remove]['score']
    # print(f"delta_c removal: {delta_c}")
    p = None
    if delta_c > 0:
        p = 1
    else:
        p = np.e ** (delta_c / t)
    step = np.random.binomial(1, p)
    if step:
        new_subgraph = cur_subslice.copy()
        new_subgraph.remove_node(leaf_to_remove)
        if len([a for a in new_subgraph.nodes if new_subgraph.nodes[a]["pertubed_node"]]) == 0:
            x=1
        return new_subgraph
    else:
        return cur_subslice.copy()


def retain_relevant_slices(G_original, module_sig_th):
    global G_modularity

    pertubed_nodes = []
    for cur_node in G_modularity.nodes():
        if G_modularity.nodes[cur_node]["pertubed_node"]:
            pertubed_nodes.append(cur_node)

    ccs = [G_modularity.subgraph(c) for c in connected_components(G_modularity)]
    params = []
    p = multiprocessing.Pool(constants.N_OF_THREADS)
    n_G_original = len(G_original)
    n_pertubed_nodes = len(pertubed_nodes)
    print(f'n_pertubed_nodes: {n_pertubed_nodes}')
    pertubed_nodes_in_ccs = []
    print(f"number of slices: {len(list(ccs))}")
    for i_cur_cc, cur_cc in enumerate(ccs):
        pertubed_nodes_in_ccs.append(
            len([cur_node for cur_node in cur_cc if G_modularity.nodes[cur_node]["pertubed_node"]]))
    perturbation_factor = min(0.7, (float(n_pertubed_nodes) / n_G_original) * (
            1 + 100 / n_G_original ** 0.5))

    for i_cur_cc, cur_cc in enumerate(ccs):
        params.append([n_G_original, cur_cc, i_cur_cc, n_pertubed_nodes, perturbation_factor])

    res = [a for a in p.map(pf_filter, params) if a is not None]
    print(f'# of slices after perturbation TH: {len(res)}/{len(params)}')
    p.close()
    if len(res) == 0:
        return nx.Graph(), [], []
    large_modules, sig_scores = zip(*res)
    fdr_bh_results = fdrcorrection0(sig_scores, alpha=module_sig_th, method='indep',
                                    is_sorted=False)

    # print(fdr_bh_results)
    # print(f'min: {min(list(fdr_bh_results[1]))}')
    passed_modules = [cur_cc for cur_cc, is_passed_th in zip(large_modules, fdr_bh_results[0]) if is_passed_th]
    return nx.algorithms.operators.union_all(passed_modules) if len(passed_modules) > 0 else nx.Graph(), [list(m.nodes)
                                                                                                          for m in
                                                                                                          passed_modules], \
           fdr_bh_results[1]

def pf_filter(params):
    global G_modularity
    n_G_original, cur_cc, i_cur_cc, n_pertubed_nodes, perturbation_factor = params
    pertubed_nodes_in_cc = [cur_node for cur_node in cur_cc if G_modularity.nodes[cur_node]["pertubed_node"]]
    if len(cur_cc) < 4 or n_pertubed_nodes == 0 or len(pertubed_nodes_in_cc) == 0:
            # or not (
            # len(pertubed_nodes_in_cc) / float(len(cur_cc)) >= perturbation_factor or len(pertubed_nodes_in_cc) / float(
            # n_pertubed_nodes) >= 0.1):
        return None
    else:
        score = hypergeom.sf(len(pertubed_nodes_in_cc), n_G_original, n_pertubed_nodes,
                             len(cur_cc)) \
                + hypergeom.pmf(len(pertubed_nodes_in_cc), n_G_original, n_pertubed_nodes,
                                len(cur_cc))
        return (cur_cc, score)

def get_final_modules(G, G_putative_modules):
    module_sigs = []
    G_putative_modules=[a for a in G_putative_modules if len(a.nodes) > 4 and len([n for n in a.nodes if G.nodes[n]["pertubed_node"]])>0]
    print(f"n of putative after filtering by size: {len(G_putative_modules)}")
    for i_cur_module, cur_G_module in enumerate(G_putative_modules):
        pertubed_nodes_in_cc = [cur_node for cur_node in cur_G_module.nodes if G.nodes[cur_node]["pertubed_node"]]
        pertubed_nodes = [cur_node for cur_node in G.nodes if G.nodes[cur_node]["pertubed_node"]]

        sig_score = hypergeom.sf(len(pertubed_nodes_in_cc), len(G.nodes), len(pertubed_nodes),
                                 len(cur_G_module.nodes)) \
                    + hypergeom.pmf(len(pertubed_nodes_in_cc), len(G.nodes), len(pertubed_nodes),
                                    len(cur_G_module.nodes))

        # final_module_threshold =  0.05 / len(G_putative_modules)
        print(sig_score)
        # if sig_score <= final_module_threshold:
        module_sigs.append((cur_G_module, sig_score )) # / len(G_putative_modules)

    module_sigs = sorted(module_sigs, key=lambda a: a[1])

    fdr_bh_results = fdrcorrection0([a[1] for a in module_sigs], alpha=0.05, method='indep',
                                    is_sorted=False)

    return sorted([a[0] for a, b in zip(module_sigs, fdr_bh_results[0]) if b],key=lambda a: -len(a))[:4]


def main(active_genes_file, network_file, slices_file=None, slice_threshold=0.3, module_threshold=0.05, prize_factor=0,
         ts=100):
    print("start running DOMINO...")
    if os.path.exists(f'{network_file}.pkl') and constants.USE_CACHE:
        G = pickle.load(open(f'{network_file}.pkl', 'rb'))
        print(f'network\' pkl is loaded: {network_file}.pkl')
    else:
        print(f'generating graph from {network_file}')
        G = build_network_static(network_file)
        pickle.dump(G, open(f'{network_file}.pkl', 'wb+'))
        print(f'network\' pkl is saved: {network_file}.pkl')

    # clean subslice from cycles and isolated nodes
    G.remove_edges_from(list(nx.selfloop_edges(G)))
    G.remove_nodes_from(list(nx.isolates(G)))

    print("done building network")
    scores = extract_scores(active_genes_file)
    G = add_scores_to_nodes(G, scores)

    modularity_connected_components = read_preprocessed_slices(slices_file)

    global G_modularity
    prune_network_by_modularity(G, modularity_connected_components, os.path.join(os.path.split(slices_file)[0],
                                                                                 os.path.split(network_file)[1].split(
                                                                                     ".")[0] + "." +
                                                                                 os.path.split(slices_file)[1].split(
                                                                                     ".")[0] + ".pkl"))
    G_modularity, relevant_slices, qvals = retain_relevant_slices(G, slice_threshold)
    print(f'{len(relevant_slices)} relevant slices were retained with threshold {slice_threshold}')
    params = []
    for i_slice, slice in enumerate(relevant_slices):
        params.append([G, G.subgraph(slice), i_slice, ts, relevant_slices, prize_factor, module_threshold])
    p = multiprocessing.Pool(constants.N_OF_THREADS)
    putative_modules = reduce(lambda a, b: a + b, p.map(sa_for_slice, params), [])
    p.close()
    print(f'n of putative modules: {len(putative_modules)}')
    final_modules = get_final_modules(G, putative_modules)
    print(f'n of final modules: {len(final_modules)} (n={[len(list(m)) for m in final_modules]})')
    return final_modules


if __name__=="__main__":
    ms=main('/home/gaga/hagailevi/omics/input/active_gene_test.txt', '/home/gaga/hagailevi/omics/input/dip.sif', slices_file='/home/gaga/hagailevi/omics/input/dip_slices.sif.pkl', slice_threshold=0.3, module_threshold=0.05,
             prize_factor=0,
             ts=100)
    print("=================")
    print([m.nodes for m in ms])