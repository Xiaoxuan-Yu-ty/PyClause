# 
import ast

from neo4j import GraphDatabase
import networkx as nx
import pickle
import pandas as pd
import re
from tqdm import tqdm
import matplotlib.pyplot as plt
import math
import argparse

import os
import sys
from utils import extract_hgnc_name, load_graph
from cmpa import compute_IF_score, run_cmpa


def extract_subgraphs_from_neo4j(driver, query=None):

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
    nodes = [{'name':dict(n)['id'],"properties": dict(n)} for n in record["nodes"]]
    
    # Extract relationships with their properties
    edges = [{
        "start": dict(r.start_node)['id'], 
        "end": dict(r.end_node)['id'], 
        "type": r.type, 
        "properties": dict(r)
    } for r in record["rels"]]
    
    return {"nodes": nodes, "edges": edges}

def get_ad_subgraph_names(subgraph_ad:dict):
    """_summary_

    Args:
        subgraph_ad (dict): _description_

    Returns:
        _type_: _description_
    """
    subgraph_names = []
    for edge_dict in subgraph_ad['edges']:
        subg = edge_dict.get('properties',{}).get('subgraph','')
        #print(type(subg))
        if isinstance(subg, str) and subg.strip():
            try:
                subg_list = ast.literal_eval(subg)
                subgraph_names.extend(subg_list)
            except:
                print(subg)
    subgraph_names = set(subgraph_names)
    print('The number of subgraphs:',len(subgraph_names))
    
    return subgraph_names

def get_ad_subgraph_query(subgraph_name):
    query = f"""
        MATCH p = (s)-[rels_in_path:increases|decreases|regulates*..4]-(o)
        WHERE ALL(n in nodes(p) WHERE n:Protein OR n:Gene) 
            AND ALL(rel in rels_in_path WHERE rel.subgraph CONTAINS '{subgraph_name}')
        
        // Use unique names for the UNWIND variables
        UNWIND nodes(p) AS individual_node
        UNWIND relationships(p) AS individual_rel
        
        WITH collect(DISTINCT individual_node) AS nodes, 
            collect(DISTINCT individual_rel) AS rels
        RETURN nodes, rels
    """
    return query

def extract_ad_subgraphs(driver, subgraph_names:set):
    """_summary_

    Args:
        driver (_type_): _description_
        subgraph_names (set): _description_

    Returns:
        _type_: _description_
    """
    subgraphs_dict = {}
    for subgraph_name in subgraph_names:
        query = get_ad_subgraph_query(subgraph_name)
        #print(query)

        # extract subgraph and store in a dict
        subg = extract_subgraphs_from_neo4j(driver, query)
        #print(f"{subgraph_name} has {len(subg['nodes'])} nodes and {len(subg['edges'])} edges")
        subg_name = '_'.join(subgraph_name.split(' '))
        subgraphs_dict[subg_name]=subg
        #break
    return subgraphs_dict

def convert_adsubgraphs_to_nx(subgraphs_dict):
    confidence_map = {'Very High':1, 'High': 0.8, 'Medium': 0.6, 'Low':0.2}
    subG_dict = {}
    for subg_name, subg in subgraphs_dict.items():
        #print(f'\nPathway name: {subg_name}')
        
        # create a subgraph for each pathway 
        sub_G = nx.MultiDiGraph()
        for node in subg['nodes']:
            #print(node)
            id = node.get('name','')
            label = node.get('properties',{}).get('type',"")
            #print(label)
            protein = node.get('properties',{}).get('name',"")
            properties = node.get('properties')
            #properties.pop('label', None)
            sub_G.add_node(id, labels=label, protein=protein, properties=properties)
        # add edges
        #print(f"Added {len(sub_G.nodes)} nodes to subGraph")
        for edge in subg['edges']:
            #print(edge)
            src = edge['start']
            dst = edge['end']
            rel = edge['type']
            confidence = edge.get('properties',{}).get('confidence', 0)
            
            props = edge.get('properties', {}).copy()
            # Remove 'confidence' from props so it's not passed twice
            props.pop('confidence', None)
            if isinstance(confidence,str):
                conf = ast.literal_eval(confidence)[0]
                confidence = confidence_map[conf]
            else:
                confidence = float(confidence)
            #print(confidence)
            sub_G.add_edge(src, dst, type=rel, confidence=confidence, **props)
            #break
        #print(f"Added {len(sub_G.edges)} edges to subGraph")

        subG_dict[subg_name] = sub_G
    print(f"Convert {len(subG_dict)} subgraph to networkx MultiDiGraphs.")
    
    return subG_dict

def plot_subgraph(G, node_label=False, edge_label=False, title="Biological Subgraph"):
    """
    Plots a NetworkX graph with color-coding for Protein/Gene nodes 
    and labels for relationship types.
    """
    plt.figure(figsize=(10, 6))
    
    # 1. Define Layout (Force-directed)
    # k regulates the distance between nodes
    pos = nx.spring_layout(G, k=0.5, seed=42)
    
    # 2. Color-code nodes based on their labels
    # We assume 'labels' was stored as a list in the node properties
    node_colors = []
    for node_id, data in G.nodes(data=True):
        label = data.get('labels', '')
        if 'Protein' in label.capitalize():
            node_colors.append('salmon')
        elif 'Gene' in label:
            node_colors.append('skyblue')
        else:
            node_colors.append('lightgrey')
    edge_colors = []
    for _,_,edtp in G.edges(data=True):
        edge_type = edtp.get('type','')
        if 'increase' in edge_type:
            edge_colors.append('black')
        elif 'decrease' in edge_type:
            edge_colors.append('blue')
        else:
            edge_colors.append('gray')
    # 3. Draw Nodes and Edges
    dynamic_size = max(10, 1000 / (len(G.nodes)**0.8))
    nx.draw_networkx_nodes(G, pos, node_size=dynamic_size, node_color=node_colors, alpha=0.9)
    nx.draw_networkx_edges(G, pos, width=1.5, alpha=0.5, edge_color=edge_colors)
    if node_label:
        # 4. Add Labels
        # Use 'name' property for nodes
        node_labels = {n:n for n in G.nodes}
        nx.draw_networkx_labels(G, pos, labels=node_labels, font_size=10, font_family="sans-serif")
    if edge_label:
        # Use 'type' property for edges (e.g., increases, decreases)
        edge_labels = {(u, v): d.get('type', '') for u, v, d in G.edges(data=True)}
        nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=8)

    plt.title(title)
    plt.axis('off')
    plt.show()

def ad_subgraph_to_rules(
                         output_dir:str,
                         expression_path:str="../Anyburl/data/adni_gene_cleaned.csv",
                         uri:str="bolt://localhost:7690", 
                         username:str="neo4j", 
                         passwaord:str="test1237", 
                         query:str|None=None,
                         G_subgraphs:dict|None = None):
                         
    driver = GraphDatabase.driver(uri, auth=(username, passwaord))
    if not query:
        query = """
        MATCH p = (s)-[:increases|decreases*..4]-(o)
        WHERE ALL(n in nodes(p) WHERE n:Protein OR n:Gene)
        UNWIND nodes(p) AS n
        UNWIND relationships(p) AS r
        WITH collect(DISTINCT n) AS nodes, collect(DISTINCT r) AS rels
        RETURN nodes, rels
        """
    if not G_subgraphs:
        # 1. Query Neo4j to Get subgraph names and then extract according subgraphs
        subgraph_flat = extract_subgraphs_from_neo4j(driver, query)
        subgraph_names = get_ad_subgraph_names(subgraph_ad=subgraph_flat)
        subgraphs_dict = extract_ad_subgraphs(driver, subgraph_names)
        # convert graph_dictionary to nx.graph
        G_subgraphs = convert_adsubgraphs_to_nx(subgraphs_dict=subgraphs_dict)

    # 2. Compute CMPA for all subgraphs
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
    save_path = os.path.join(output_dir, 'ad_subgraph_features.csv')
    df_patient_cmpa.to_csv(save_path, index=False)
    return df_patient_cmpa

def main():
    parser = argparse.ArgumentParser(description="Extract subgraphs from Neo4j and convert them to features by CMPA")
    parser.add_argument("--expression_path", type=str, default="../Anyburl/data/adni_gene_cleaned.csv")
    parser.add_argument("--uri", type=str, default="bolt://localhost:7690")
    parser.add_argument("--username", type=str, default="neo4j")
    parser.add_argument("--password", type=str, default="test1237")
    parser.add_argument("--output_dir", type=str, default="../results")

    args = parser.parse_args()

    df_features = ad_subgraph_to_rules(output_dir=args.output_dir,
                                       expression_path=args.expression_path,
                                       uri=args.uri,
                                       username=args.username,
                                       passwaord=args.password,
                                       query=None
                                    )



if __name__ == "__main__":
    main()