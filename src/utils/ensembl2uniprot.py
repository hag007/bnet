import os
from src import constants


u2e_dict = None
e2u_dict = None
ensembl2uniprot_dict = None
uniprot2ensembl_dict = None

def load_gene_dictionary(gene_list_file_name, gene_list_path=None):
    if gene_list_path == None:
        gene_list_path = os.path.join(constants.dir_path,"data",gene_list_file_name)
    f = open(gene_list_path,'r')
    lines = [l.strip() for l in f]
    f.close()
    return lines


def get_uniprot2ensembl_dictionary():
    lines_dict = load_gene_dictionary(constants.ENSEMBL_TO_UNIPROT)

    uniprot2ensembl = {}
    for cur in lines_dict:
        splited_line = cur.split()
        if len(splited_line) != 2: continue
        if splited_line[0].find('.') > 0:
            limit = splited_line[0].find('.')
        else:
            limit = len(splited_line[0])
            uniprot2ensembl[splited_line[1]] = splited_line[0][:limit]
    return uniprot2ensembl

def get_ensembl2uniprot_dictionary():

    lines_dict = load_gene_dictionary(constants.ENSEMBL_TO_UNIPROT)
    ensembl2gene_symbols = {}
    for cur in lines_dict:
        splited_line = cur.split()
        if len(splited_line) !=2: continue
        if splited_line[0].find('.') > 0:
            limit = splited_line[0].find('.')
        else:
            limit = len(splited_line[0])
        ensembl2gene_symbols[splited_line[0][:limit]] = splited_line[1]
    return ensembl2gene_symbols


def ensembl2uniprot_convertor(e_ids):
    global ensembl2uniprot_dict
    if ensembl2uniprot_dict is None:
        ensembl2uniprot_dict = get_ensembl2uniprot_dictionary()
    results = []
    for cur in e_ids:
        if cur.split(".")[0] in ensembl2uniprot_dict:
            results.append(ensembl2uniprot_dict[cur.split(".")[0]])
    return results


def uniprot2ensembl_convertor(uniprot_ids):
    global uniprot2ensembl_dict
    if uniprot2ensembl_dict is None:
        uniprot2ensembl_dict = get_uniprot2ensembl_dictionary()
    results = []
    for cur in uniprot_ids:
        cur=str(cur)
        if uniprot2ensembl_dict.has_key(cur):
            results.append(uniprot2ensembl_dict[cur])
    return results
