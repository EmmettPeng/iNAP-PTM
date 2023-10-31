# Obtain potentially transferable metabolites

from libsbml import readSBML
import networkx as nx
import pandas as pd
from xml.dom import minidom

def buildDG(sbml):
    '''
    Usage: reads SBML file, parses reaction and product list
    Returns: networkx directed graph
    '''

    # initate empty directed graph
    DG = nx.DiGraph()

    document = readSBML(sbml)
    model = document.getModel()

    # get reactions
    rxn = (model.getListOfReactions())

    # get reaction and product and construct directed graph
    for r in rxn:
        react = [i.getSpecies() for i in r.getListOfReactants()]
        prod = [j.getSpecies() for j in r.getListOfProducts()]
        for r in react:
            for p in prod:
                DG.add_edge(r,p)

    return DG

def getSeedSet(DG, maxComponentSize = 5):
    '''
    Usage: takes input networkX directed graph
    Returns: SeedSet dictionary{seedset:confidence score}
    Implementation follows literature description, 
    Improves upon NetCooperate module implementation which erroneously discards certian cases of SCCs (where a smaller potential SCC lies within a larger SCC)
    '''

    # get SCC
    SCC = nx.strongly_connected_components(DG)

    SeedSetConfidence = dict()

    for cc in SCC:
        # convert set to list
        cc_temp = list(cc)

        # filter out CC larger than threshold
        if len(cc_temp) > int(maxComponentSize):
            continue

        # check single element SCC
        elif len(cc_temp) == 1:
            if DG.in_degree(cc_temp[0]) == 0:
                SeedSetConfidence[cc_temp[0]] = 1.0    

        # check 2 to max threshold SCC
        else:
            #check if no out nodes
            for node in cc_temp:
                # check every edge of SCC
                for edge in DG.in_edges(node):
                    # if SCC is not self contained, then it is not considered seed set
                    if edge[0] not in cc_temp:
                        cc_temp = []
            for node in cc_temp:
                SeedSetConfidence[node] = 1/len(cc_temp) 
    
    SeedSet = SeedSetConfidence.keys()

    nonSeedSet = list(set(DG.nodes()) - set(SeedSet))

    return(SeedSetConfidence, SeedSet, nonSeedSet) 

def get_ptm(A_conf, B_nonseedset):
    '''
    Usage: Get potentially transferable metabolites from donor to receptor; A: Receptor, B: Donor
    Returns: [ptm_bigg_ids, ...] (list)
    '''

    SeedA = set(A_conf.keys())
    #SeedB = set(B_conf.keys())
    nonSeedB = set(B_nonseedset)
    #SetB = (SeedB|nonSeedB)
    # get intersects intersect (A n nonB) &! B
    intersect_seedA_nonseedB = SeedA.intersection(nonSeedB)
    #intersect_seedA_nonseedB : ptms

    return list(intersect_seedA_nonseedB)

def main():
    A_filepath, B_filepath = "MAG258.xml","MAG354.xml"
    A_graph, B_graph = buildDG(A_filepath), buildDG(B_filepath)
    A_model_name = minidom.parse(A_filepath).documentElement.getElementsByTagName("model")[0].getAttribute("metaid")
    B_model_name = minidom.parse(B_filepath).documentElement.getElementsByTagName("model")[0].getAttribute("metaid")
    A_conf, A_seedset, A_nonseedset = getSeedSet(A_graph, 5)
    B_conf, B_seedset, B_nonseedset = getSeedSet(B_graph, 5)

    # A as donor and B as receptor
    ptm_a_to_b = get_ptm(B_conf, A_nonseedset)
    ptm_b_to_a = get_ptm(A_conf, B_nonseedset)

    # Key: bigg_id
    # Value0: universal_bigg_id; Value1: name; Value2: model_list; Value3: database_links; Value4: old_bigg_ids 
    bigg_db = pd.read_csv('bigg_models_metabolites.txt', sep = '\t')
    bigg_dict = pd.Series(bigg_db.set_index('bigg_id').T.to_dict('list'))

    with open("test_output.txt", "w") as out:
        out.write('Donor\tReceptor\tbigg_id\tuniversal_bigg_id\tmetabolite_name\tbigg_model_list\tdatabase_links\n')
        for name in ptm_a_to_b:
            name = name[2:]
            out.write(A_model_name+'\t'+B_model_name+'\t'+name+'\t'+str(bigg_dict[name][0])+'\t'+str(bigg_dict[name][1])+'\t'+str(bigg_dict[name][2])+'\t'+str(bigg_dict[name][3])+'\n')
        for name in ptm_b_to_a:
            name = name[2:]
            out.write(B_model_name+'\t'+A_model_name+'\t'+name+'\t'+str(bigg_dict[name][0])+'\t'+str(bigg_dict[name][1])+'\t'+str(bigg_dict[name][2])+'\t'+str(bigg_dict[name][3])+'\n')

if __name__ == '__main__':
    main()