import argparse

from neo4j import GraphDatabase
import pandas as pd
from tqdm import tqdm
import os 
import sys
try:
    base_dir = os.path.dirname(os.path.abspath(__file__))
except NameError:
    base_dir = os.getcwd()
print(f'base_dir is {base_dir}')
sys.path.append(os.path.dirname(base_dir))
from util.cmpa import (
                       run_cmpa, 
                       extract_subgraphs_from_neo4j, 
                       convert_graphdict_to_nx,
                       subgraph_to_rules)
import networkx as nx



def extract_bp_subgraphs_dedup(G, hop=2, min_size=5):
    bp_nodes = [n for n, d in G.nodes(data=True) if d.get("labels") == "BiologicalProcess"]
    
    subgraphs = {}
    protein_coverage = {}  # track which BPs cover each protein
    
    for bp in bp_nodes:
        neighbors = nx.single_source_shortest_path_length(G, bp, cutoff=hop)
        sg = G.subgraph(list(neighbors.keys())).copy()
        
        if sg.number_of_nodes() < min_size:
            continue
        
        subgraphs[bp] = sg
        
        # Track protein membership
        for node in sg.nodes():
            if G.nodes[node].get("labels") == "Protein":
                protein_coverage.setdefault(node, []).append(bp)
    # Report bp-subgrphs stats
    subgraph_num_nodes = [subg.number_of_nodes() for subg in subgraphs.values()]
    print(f"Number of subgraphs: {len(subgraphs)}")
    print(f"Subgraph node count: min={min(subgraph_num_nodes)}, max={max(subgraph_num_nodes)}")
    # Report overlap stats
    overlap = {p: bps for p, bps in protein_coverage.items() if len(bps) > 1}
    print(f"\nProteins appearing in multiple BP subgraphs: {len(overlap)}")
    print(f"Max overlap: {max(len(v) for v in protein_coverage.values())} BPs")
    
    return subgraphs, protein_coverage

def subgraph_to_rules(
                    G_subgraphs:dict,
                    output_dir:str,
                    expression_path:str="../Anyburl/data/adni_gene_cleaned.csv",
                    ):
    # 1. Compute CMPA for all subgraphs
    exp_df = pd.read_csv(expression_path, index_col=0).T
    df_patient_cmpa = pd.DataFrame(index=exp_df.index, columns = list(G_subgraphs.keys()))
    
    for i in tqdm(range(len(exp_df)), desc="Computing Subgraph CMPA Score For Patients"):
        patient_data = exp_df.iloc[i].to_dict()

        df_subgraphs = {}
        subgraph_scores = []
        for subg_name, subG in G_subgraphs.items():
            #plot_subgraph(subG, True, True, subg_name)
            print('='*60)
            print(f"Run CMPA on {subg_name}")
            df_subG = run_cmpa(patient_data, subG)
            df_subgraphs[subg_name] = df_subG
            
            subgraph_scores.append(float(df_subG[df_subG['hub'] == subg_name]['score']))
        df_patient_cmpa.iloc[i] = subgraph_scores
    
    # save features file
    #save_path = os.path.join(output_dir, 'ad_subgraph_features.csv')
    df_patient_cmpa.to_csv(output_dir)
    return df_patient_cmpa


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--KG", type=str, default='health',
                        choices=['ad', 'health'])
    args = parser.parse_args()
    # Connection to Neo4j
    driver_ad = GraphDatabase.driver("bolt://localhost:7690", auth=('neo4j','test1237'))
    driver_health = GraphDatabase.driver("bolt://localhost:7688", auth=('neo4j','test12345'))

    query_BP = """
            MATCH p=(s)-[]->(o) 
            UNWIND nodes(p) AS n
            UNWIND relationships(p) AS r
            WITH collect(DISTINCT n) AS nodes, collect(DISTINCT r) AS rels
            RETURN nodes, rels
            """
    kg = args.KG
    if kg == 'ad':
        # convert Neo4j cypher to networkx graph
        BP_dict = extract_subgraphs_from_neo4j(driver=driver_ad, 
                                                    query=query_BP, 
                                                    ad=True)
        bp_G = convert_graphdict_to_nx(graph_dict=BP_dict, ad=True)
        
        # remove nodes that are not Proteins or BiologicalProcess
        node_to_remove = [n for n, attr in bp_G.nodes(data=True) if attr.get('labels') not in ['Protein', 'BiologicalProcess']]
        len(node_to_remove)
        bp_G.remove_nodes_from(node_to_remove)
        print(f"nodes:{bp_G.number_of_nodes()}, edges:{bp_G.number_of_edges()}")

        # get BP-Ego graphs
        subgraphs, protein_coverage = extract_bp_subgraphs_dedup(G=bp_G,
                                                                hop=2,
                                                                min_size=5)
        df_bp_cmpa = subgraph_to_rules(G_subgraphs=subgraphs,
                            output_dir = "../results/BPsubgraphs_ad.csv",
                            expression_path = "../Anyburl/data/adni_gene_cleaned.csv")
    elif kg=='health':
        # convert Neo4j cypher to networkx graph
        BP_dict = extract_subgraphs_from_neo4j(driver=driver_health, 
                                                    query=query_BP, 
                                                    ad=False)
        bp_G = convert_graphdict_to_nx(graph_dict=BP_dict, ad=False)

        # change label
        for node, attr in bp_G.nodes(data=True):
            if node.startswith('bp'):
                label = attr['labels']
                new_label = ''.join([word.capitalize() for word in label.split('_')])
                #print(new_label)
                attr['labels'] = new_label
        
        # remove nodes that are not Proteins or BiologicalProcess
        node_to_remove = [n for n, attr in bp_G.nodes(data=True) if attr.get('labels') not in ['Protein', 'BiologicalProcess']]
        len(node_to_remove)
        bp_G.remove_nodes_from(node_to_remove)
        print(f"nodes:{bp_G.number_of_nodes()}, edges:{bp_G.number_of_edges()}")

        # get BP-Ego graphs
        subgraphs, protein_coverage = extract_bp_subgraphs_dedup(G=bp_G,
                                                                hop=2,
                                                                min_size=5)
        df_bp_cmpa = subgraph_to_rules(G_subgraphs=subgraphs,
                            output_dir = "../results/BPsubgraphs_health.csv",
                            expression_path = "../Anyburl/data/adni_gene_cleaned.csv")



if __name__ == "__main__":
    main()
    