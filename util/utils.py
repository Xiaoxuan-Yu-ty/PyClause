
import re

from neo4j import GraphDatabase
import networkx as nx
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle

def load_graph(kg_path)-> nx.MultiDiGraph:
    """Load nx.graph from .pkl file

    Args:
        kg_path (str): _description_

    Returns:
        nx.MultiDigraph: _description_
    """
    with open(kg_path, 'rb') as f:
        kg = pickle.load(f)
    print(f"Graph has {kg.number_of_nodes()} nodes and {kg.number_of_edges()} edges.")
    return kg

def get_nx_graph_from_neo4j(output:str,
                            url:str="bolt://localhost:7688", 
                            user:str="neo4j", 
                            passwd:str="test12345")->nx.MultiDiGraph:
    """Extract KG from Neo4j and convert it to networkx.MultiDiGraph

    Args:
        output (str): filename to save the graph in .pkl format
        url (str, optional): url to the Neo4j browser. Defaults to "bolt://localhost:7688".
        user (str, optional): user name. Defaults to "neo4j".
        passwd (str, optional): password. Defaults to "test12345".

    Returns:
        nx.MultiDiGraph: the converted graph
    """
    driver = GraphDatabase.driver(url, auth=(user, passwd))
    G = nx.MultiDiGraph() # MultiDiGraph handles multiple edges between nodes
    with driver.session() as session:
        # Match everything
        result = session.run("MATCH (n)-[r]->(m) RETURN n, r, m")
        for record in result:
            u = record['n'].element_id
            v = record['m'].element_id
            # Add nodes and edges with their properties
            G.add_node(u, **dict(record['n']))
            G.add_node(v, **dict(record['m']))
            G.add_edge(u, v, key=record['r'].element_id, type=record['r'].type, **dict(record['r']))
    
    with open(output, 'wb') as f:
        pickle.dump(G, f)
    print(f"Successfully created NetworkX graph with {len(G.nodes)} nodes!")
    return G

def add_node(tx, node_id, attrs):
    """Write nodes to Neo4j Database"""
    node_label = attrs.get('label', 'Entity')
    query = (
        f"MERGE (n:`{node_label}` {{id: $id}}) "
        "SET n += $props"
    )
    
    # n += $props sets all keys in the dict as properties at once
    tx.run(query, id=str(node_id), props=attrs)

def add_relationship(tx, head, tail, rel_type, attrs):
    """Write Relations to Neo4j Database"""
    
    # Match the existing nodes by ID
    query = (
        "MATCH (a {id: $head}), (b {id: $tail}) "
        f"CREATE (a)-[r:`{rel_type}`]->(b) "
        "SET r.confidence = $confidence, "
        "    r.subgraph = $subgraph, "
        "    r.evidence = $evidence,"
        "    r.disease = $disease,"
        "    r.pmid = $pmid"
    )
    
    tx.run(query, 
           head=str(head), 
           tail=str(tail), 
           confidence=attrs.get('annotation', {}).get('confidence'),
           subgraph=attrs.get('annotation', {}).get('subgraph'),
           evidence=attrs.get('evidence'),
           disease=attrs.get('disease'),
           pmid = attrs.get('pmid'))

def add_kg_to_neo4j(G:nx.MultiDiGraph,
                        url:str="bolt://localhost:7689", 
                        user:str="neo4j", 
                        passwd:str="test1236"):
    """Write KG Triples to Neo4j Database.

    Args:
        G (nx.MultiDiGraph): _description_
        url (_type_, optional): _description_. Defaults to "bolt://localhost:7689".
        user (str, optional): _description_. Defaults to "neo4j".
        passwd (str, optional): _description_. Defaults to "test1236".
    """
    driver = GraphDatabase.driver(uri=url, auth=(user, passwd))

    with driver.session() as session:
        # add nodes:
        for node_id, attrs in G.nodes(data=True):
            session.execute_write(add_node, node_id, attrs)
        # add edges
        for head, tail, rel_type, attrs in G.edges(data=True, keys=True):
            session.execute_write(add_relationship, head, tail, rel_type, attrs)
    # close driver
    driver.close()


def extract_uniprotkb_name_species(text:str, pattern:str = r'(UniProtKB:"(\w+)_(\w+)")'):
    match = re.search(pattern, text)
    if match:
        #print(match.group())
        name = match.group(2)
        species = match.group(3)
        return name, species
    else:
        print('No match found.')

def extract_hgnc_name(text, pattern:str=r'(HGNC:"(\w+)")'):
    match = re.search(pattern, text)
    if match:
        name = match.group(2)
        #print(name)
        return name
    else:
        return None
