import math
import random
import sys

sys.path.insert(0, '../')
import os
import numpy as np

import pandas as pd

from fastsemsim.semsim import *
import fastsemsim.ontology as  ontologies
from fastsemsim.ac.AnnotationCorpus import AnnotationCorpus
from fastsemsim.semsim.SetSemSim import SetSemSim
import fastsemsim

import constants

import multiprocessing

from utils.daemon_multiprocessing import func_star
from src.utils.ensembl2uniprot import ensembl2uniprot_convertor
import simplejson as json
from simplejson.errors import JSONDecodeError

has_go_resources = False
has_go_metadata = False
has_go_hierarchy = False
loading_metadata = False
ROOT = 'GO:0008150'
CUTOFF = 7
ssbatch= None

def init_go_metadata():
    global has_go_metadata
    global loading_metadata
    global ssbatch
    global ontology
    global ac

    if not has_go_metadata and not loading_metadata:
        loading_metadata = True
        ontology_type = 'GeneOntology'
        ignore_parameters = {'ignore': {}}
        source_type = 'obo'
        source = os.path.join(os.path.join(constants.GO_DIR, constants.GO_FILE_NAME))

        print("Loading GO files..")

        ontology = ontologies.load(source_file=source, ontology_type=ontology_type,
                                   parameters=ignore_parameters)

        ac = AnnotationCorpus(ontology)
        ac.parse(os.path.join(constants.GO_DIR, "goa_human.gaf"), "gaf-2.0")
        ac.isConsistent()

        has_go_metadata = True
        loading_metadata = False

        ssbatch = fastsemsim.init_batchsemsim(ontology=ontology, ac=ac, semsim_type='obj', semsim_measure='Resnik',
                                              mixing_strategy='BMA')
        ssbatch.set_root('biological_process')
        ssbatch.set_output(output=None)


def init_go_hierarchy():
    global has_go_hierarchy

    if not has_go_hierarchy:
        global vertices

        import utils.go_hierarchy as go_hierarchy
        ROOT = 'GO:0008150'
        dict_result, go2geneids, geneids2go, entrez2ensembl = go_hierarchy.build_hierarchy(
            roots=[ROOT])
        vertices = list(dict_result.values())[0]['vertices']

        has_go_hierarchy = True


def init_go_resources():
    global has_go_resources

    if not has_go_resources:

        if not has_go_metadata:
            init_go_metadata()

        if not has_go_hierarchy:
            init_go_hierarchy()

        has_go_resources = True


semsim_dict = {}


def get_semsim(sim_method):
    if sim_method not in semsim_dict:
        init_go_metadata()
        semsim_dict[sim_method] = SetSemSim(ontology, ac, TSS=sim_method, MSS="BMA")

    return semsim_dict[sim_method]


def calc_similarity(mat_adj, x, y, semsim):
    key = "{}_{}".format(x, y)
    key_inv = "{}_{}".format(y, x)
    if mat_adj[key] != -200: return
    if (x == y):
        mat_adj[key] = -100
        return

    mat_adj[key] = get_semsim(semsim).SemSim(x, y)
    # if mat_adj[key]<0:
    # print i_x,i_y,mat_adj[key]

    # print(mat_adj[key])
    if mat_adj[key] is None or np.isnan(mat_adj[key]):
        mat_adj[key] = -100
    mat_adj[key_inv] = mat_adj[key]


def calc_similarity_pairs(set_0, set_1):
    sm=0
    for cur_0 in set_0:
        for cur_1 in set_1:
            if len(ensembl2uniprot_convertor([cur_0]) + ensembl2uniprot_convertor([cur_1])) < 2:
                sm+=0.5
            else:
                cur=ssbatch.SemSim(query=[ensembl2uniprot_convertor([cur_0]) + ensembl2uniprot_convertor([cur_1])]).iloc[0,2]
                if cur is None:
                    sm +=0.5
                else:
                    sm+=cur
    return 1/(1+(sm/(len(set_0)*len(set_1))))

def calc_similarity_matrix(set_0, set_1, pf, cache_file, sim_method='Resnik'):
    cache_loaded = False
    if os.path.exists(cache_file):
        try:
            adj = json.load(open(cache_file, 'r'))
            cache_loaded = True
        except:
            pass

    if not cache_loaded:

        init_go_resources()

        semsim = get_semsim(sim_method)
        manager = multiprocessing.Manager()
        adj = manager.dict()
        for x in set_0:
            for y in set_1:
                adj["{}_{}".format(x, y)] = -200
        params = []
        for i_x, x in enumerate(set_0):
            for i_y, y in enumerate(set_1):
                calc_similarity(adj, x, y, semsim)
                # params.append([calc_similarity, [adj, i_x, i_y, x, y, semsim]])
        print("len(params): {}".format(len(params)))

        # p = multiprocessing.Pool(pf)
        # p.map(func_star, params)
        # p.close()
        # p.join()
        for p in params:
            p[0](*p[1])

    open(cache_file, 'w+').write(json.dumps(dict(adj)))
    return adj


def calc_intra_similarity(all_go_terms, pf, enrichment_scores, cache_file, semsim, cutoff=CUTOFF, set_0=[], set_1=[],
                          reduce_list=True):
    if not has_go_hierarchy:
        init_go_hierarchy()

    cache_loaded = True
    try:
        adj = json.load(open(cache_file, 'r'))
    except Exception:
        cache_loaded = False

    if cache_loaded:
        if all_go_terms is None:
            all_go_terms_r = list(np.unique(np.array([cur.split("_") for cur in adj]).flatten()))
        else:
            all_go_terms_r = list(all_go_terms)

    else:
        all_go_terms_r = list(all_go_terms)
        manager = multiprocessing.Manager()
        adj = manager.dict()
        for x in all_go_terms_r:
            for y in all_go_terms_r:
                adj["{}_{}".format(x, y)] = -200
        params = []
        for i_x, x in enumerate(all_go_terms_r):
            for i_y, y in enumerate(all_go_terms_r):
                # print x, y
                # calc_similarity(adj, x, y, semsim)
                params.append([calc_similarity, [adj, x, y, semsim]])
        print("len(params): {}".format(len(params)))

        init_go_metadata()

        #         p = multiprocessing.Pool(pf)
        #         p.map(func_star, params)
        #         p.close()
        #         p.join()
        for p in params:
            p[0](*p[1])
        adj = dict(adj)
        json.dump(adj, open(cache_file, 'w+'))

    all_go_terms_o = list(all_go_terms_r)

    if reduce_list:
        execute_revigo(adj, all_go_terms_o, all_go_terms_r, enrichment_scores, cutoff, set_0, set_1)

    # print(all_go_terms_r)
    return all_go_terms_r, all_go_terms_o, adj


def execute_revigo(adj, all_go_terms_o, all_go_terms_r, enrichment_scores, cutoff=CUTOFF, set_1=[], set_2=[]):
    for cur_go_1, cur_go_id_1 in enumerate(all_go_terms_o):
        for cur_go_2, cur_go_id_2 in enumerate(all_go_terms_o):

            if not "{}_{}".format(cur_go_id_1, cur_go_id_2) in adj: continue

            if cur_go_id_1 not in all_go_terms_r or cur_go_id_2 not in all_go_terms_r or cur_go_id_2 == cur_go_id_1:
                continue
            if adj["{}_{}".format(cur_go_id_1, cur_go_id_2)] > cutoff:

                if cur_go_id_1 in set_1 and len(set_1) < len(set_2):
                    all_go_terms_r.remove(cur_go_id_1)
                elif cur_go_id_1 in set_2 and len(set_2) < len(set_1):
                    all_go_terms_r.remove(cur_go_id_1)
                elif cur_go_id_2 in set_1 and len(set_1) < len(set_2):
                    all_go_terms_r.remove(cur_go_id_2)
                elif cur_go_id_2 in set_2 and len(set_2) < len(set_1):
                    all_go_terms_r.remove(cur_go_id_2)

                elif enrichment_scores.loc[cur_go_id_1, 'hg_pval_max'] > enrichment_scores.loc[
                    cur_go_id_2, 'hg_pval_max'] + 1:
                    all_go_terms_r.remove(cur_go_id_2)
                elif enrichment_scores.loc[cur_go_id_1, 'hg_pval_max'] + 1 < enrichment_scores.loc[
                    cur_go_id_2, 'hg_pval_max']:
                    all_go_terms_r.remove(cur_go_id_1)
                    break
                elif is_parent(cur_go_id_1, cur_go_id_2):
                    all_go_terms_r.remove(cur_go_id_2)
                elif is_parent(cur_go_id_2, cur_go_id_1):
                    all_go_terms_r.remove(cur_go_id_1)
                else:
                    all_go_terms_r.remove([cur_go_id_1, cur_go_id_2][int(math.floor(random.random() * 2))])


def is_parent(parent, child):
    c_parents = vertices[child]["obj"].parents
    if child == ROOT:
        return False

    elif len(c_parents) == 0:
        return True

    else:
        return any([is_parent(parent, a.id) for a in list(c_parents)])

    return adj, all_go_terms


init_go_metadata()