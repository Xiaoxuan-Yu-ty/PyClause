# 
from neo4j import GraphDatabase
import networkx as nx
import pickle
import pandas as pd
import re
from tqdm import tqdm
import matplotlib.pyplot as plt
import math
import ast
import igraph as ig
import leidenalg

import os
import sys
from cmpa import compute_IF_score, run_cmpa

# Connection to Neo4j
driver_ad = GraphDatabase.driver("bolt://localhost:7690", auth=('neo4j','test1237'))
driver_health = GraphDatabase.driver("bolt://localhost:7688", auth=('neo4j','test12345'))

def extract_subgraphs_from_neo4j(driver, query=None, ad=False):
    if not query:
        query = """
        MATCH p = (s)-[:increases|decreases*..4]-(o)
        WHERE ALL(n in nodes(p) WHERE n:Protein OR n:Gene)
        UNWIND nodes(p) AS n
        UNWIND relationships(p) AS r
        WITH collect(DISTINCT n) AS nodes, collect(DISTINCT r) AS rels
        RETURN nodes, rels
        """
    def query_neo4j(tx):
        result = tx.run(query)
        return result.single()
    with driver.session() as session:
        record = session.execute_read(query_neo4j)
    
    # Extract nodes with their properties
    if not ad:
        nodes = [{'name':dict(n)['bel'],"properties": dict(n)} for n in record["nodes"]]
        # Extract relationships with their properties
        edges = [{
            "start": dict(r.start_node)['bel'], 
            "end": dict(r.end_node)['bel'], 
            "type": r.type, 
            "properties": dict(r)
        } for r in record["rels"]]
    else:
        nodes = [{'name':dict(n)['id'],"properties": dict(n)} for n in record["nodes"]]
        # Extract relationships with their properties
        edges = [{
            "start": dict(r.start_node)['id'], 
            "end": dict(r.end_node)['id'], 
            "type": r.type, 
            "properties": dict(r)
        } for r in record["rels"]]
    
    return {"nodes": nodes, "edges": edges}


def convert_graphdict_to_nx(graph_dict, ad=False):
    
    # can serve as edge_weight
    confidence_map = {'Very High':1, 'High': 0.8, 'Medium': 0.6, 'Low':0.3}
      
    G = nx.MultiDiGraph()
    
    # add nodes to G
    for node in graph_dict['nodes']:
        #print(node)
        id = node.get('name','')
        label = node.get('properties',{}).get('type',"")
        #print(label)
        protein = node.get('properties',{}).get('name',"")
        protein = protein.split('_')[0]
        properties = node.get('properties')
        #properties.pop('label', None)
        G.add_node(id, labels=label, protein=protein, properties=properties)
    print(f"Added {len(G.nodes)} nodes to subGraph")
        
    # add edges
    for edge in graph_dict['edges']:
        #print(edge)
        src = edge['start']
        dst = edge['end']
        rel = edge['type']
        props = edge.get('properties', {}).copy()
        if ad:
            confidence = edge.get('properties',{}).get('confidence', 0.5)
            #print(confidence)
            # Remove 'confidence' from props so it's not passed twice
            props.pop('confidence', None)
            if isinstance(confidence,str):
                conf = ast.literal_eval(confidence)[0]
                confidence = confidence_map[conf]
            elif isinstance(confidence, list):
                conf = ast.literal_eval(confidence[0])[0]
                if isinstance(conf, str):
                    confidence = confidence_map[conf]
            else:
                confidence = float(confidence)
            #print(confidence)
        else:
            confidence = 1.0 if props.get('pmid',0) != 0 else 0.0 
        G.add_edge(src, dst, type=rel, weight=confidence, **props)
        #break
    print(f"Added {len(G.edges)} edges to subGraph")

    return G


def networkx_to_igraph(G):
    """Convert weighted NetworkX graph to igraph."""
    nodes = list(G.nodes())
    node_idx = {n: i for i, n in enumerate(nodes)}
    
    edges = [(node_idx[u], node_idx[v]) for u, v in G.edges()]
    weights = [d.get("weight", 1.0) for _, _, d in G.edges(data=True)]
    
    ig_graph = ig.Graph(n=len(nodes), edges=edges)
    ig_graph.vs["name"] = nodes
    ig_graph.es["weight"] = weights
    return ig_graph, nodes

def leiden_partition(G, resolution=1.0):
    """
    resolution > 1.0 → smaller, more communities
    resolution < 1.0 → larger, fewer communities
    """
    ig_G, nodes = networkx_to_igraph(G)
    partition = leidenalg.find_partition(
        ig_G,
        leidenalg.RBConfigurationVertexPartition,  # supports weights + resolution
        weights="weight",
        resolution_parameter=resolution,
        n_iterations=10,
        seed=42
    )
    
    # Map back to original node names
    communities = []
    for community in partition:
        communities.append([nodes[i] for i in community])
    
    return communities

def extract_feature_subgraphs(G, communities, min_size=5, max_size=70):
    """
    Extract subgraphs from communities, splitting oversized ones recursively.
    """
    feature_subgraphs = []
    
    for comm_nodes in communities:
        subgraph = G.subgraph(comm_nodes).copy()
        
        if len(comm_nodes) < min_size:
            pass  # too small, discard
        
        elif len(comm_nodes) <= max_size:
            feature_subgraphs.append(subgraph)
        
        else:
            # Recursively re-partition oversized communities
            sub_communities = leiden_partition(subgraph, resolution=2.0)
            for sub_comm in sub_communities:
                if len(sub_comm) >= min_size:
                    feature_subgraphs.append(G.subgraph(sub_comm).copy())
    
    return feature_subgraphs

def neo4j_to_featureGraphsDict(driver, G=None, query=None, max_size=70, min_size=5):

    if not G:
        # 1. export neo4j query graph to networkx
        graph_dict = extract_subgraphs_from_neo4j(driver,query=query)
        G = convert_graphdict_to_nx(graph_dict)

    # 2. Find all connected components, sorted by size
    components = sorted(nx.connected_components(G.to_undirected()), key=len, reverse=True)
    print(f"Number of components: {len(components)}")

    subgraphs_dict = {}
    for i, c in enumerate(components):
        if len(c) < min_size:
            subgraphs_dict[f'subgraph_{i}'] = G.subgraph(components[i]).copy()
        #print(f"  Component {i}: {len(c)} nodes")

    # Extract the giant component as its own subgraph
    giant = G.subgraph(components[0]).copy()
    print(f"\nGiant component: {giant.number_of_nodes()} nodes, {giant.number_of_edges()} edges")

    # 3. Paritition Giant componnets to smaller ones
    communities = leiden_partition(giant, resolution=1.0)
    # 4. filter partitioned Giant-subgraph and store subgraphs to dict
    feature_subgraphs = extract_feature_subgraphs(giant, communities, min_size=5, max_size=80)
    for i, sg in enumerate(feature_subgraphs):
        subgraphs_dict[f'subgraph_{i+len(components)}'] = sg
        #print(f"  Subgraph {i}: {sg.number_of_nodes()} nodes, {sg.number_of_edges()} edges")

    print(f"\nFinal feature subgraphs: {len(subgraphs_dict)}")

    return subgraphs_dict

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
            
            subgraph_scores.append(df_subG['score'].sum())
        df_patient_cmpa.iloc[i] = subgraph_scores
    
    # save features file
    #save_path = os.path.join(output_dir, 'ad_subgraph_features.csv')
    df_patient_cmpa.to_csv(output_dir)
    return df_patient_cmpa

def main():
    # Connection to Neo4j
    driver_ad = GraphDatabase.driver("bolt://localhost:7690", auth=('neo4j','test1237'))
    driver_health = GraphDatabase.driver("bolt://localhost:7688", auth=('neo4j','test12345'))

    query_ad = """
        MATCH p = (s)-[:increases|decreases*..4]-(o)
        WHERE ALL(n in nodes(p) WHERE n:Protein OR n:Gene)
        UNWIND nodes(p) AS n
        UNWIND relationships(p) AS r
        WITH collect(DISTINCT n) AS nodes, collect(DISTINCT r) AS rels
        RETURN nodes, rels
        """
    query_health = """
        MATCH p = (s)-[:increases|decreases*..4]-(o)
        WHERE ALL(n in nodes(p) WHERE n:protein OR n:gene)
        UNWIND nodes(p) AS n
        UNWIND relationships(p) AS r
        WITH collect(DISTINCT n) AS nodes, collect(DISTINCT r) AS rels
        RETURN nodes, rels
        """
    query = """
        "MATCH p=(s)-[]->(o) 
        UNWIND nodes(p) AS n
        UNWIND relationships(p) AS r
        WITH collect(DISTINCT n) AS nodes, collect(DISTINCT r) AS rels
        RETURN nodes, rels
        """
    feature_graphs_dict = neo4j_to_featureGraphsDict(driver_ad, query=query)
    df_ad_cmpa = subgraph_to_rules(feature_graphs_dict,'../results/ad_clutser_features.csv')

    feature_graphs_dict = neo4j_to_featureGraphsDict(driver_health, query=query_health)
    df_health_cmpa = subgraph_to_rules(feature_graphs_dict,'../results/health_clutser_features.csv')

if __name__ == "__main__":
    main()