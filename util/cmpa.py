# CMPA Algorithm - Cyclic Graph Compatible
import ast

import pandas as pd
import networkx as nx
from neo4j import GraphDatabase
import matplotlib.pyplot as plt
import math

from tqdm import tqdm


def compute_IF_score(graph, node):
    IF_score = graph.nodes[node]['base_score'] # always startwith base score instead of IF_score, otherwise not converge
    
    for pred in graph.predecessors(node):
        edge_dict = graph[pred][node]
        #print(edge_dict)
        rel_type = graph[pred][node][0]['type']
        if 'increases' in rel_type.lower():
            # print(pred, "-->", node)
            # print("Before: ", graph.nodes[node]["IF_score"])
            
            IF_score += graph.nodes[pred]["IF_score"] * 1.0
            
            # print("After: ", IF_score)
        elif 'decreases' in rel_type.lower():
            # print(pred, "--|", node)
            # print("Before: ", graph.nodes[node]["IF_score"])
            
            IF_score += graph.nodes[pred]["IF_score"] * -1.0
            
            # print("After: ", IF_score)
        else:
            IF_score += graph.nodes[pred]["IF_score"] * 0.1
    
    return IF_score

def get_max_iter(graph):
    # Diameter of the underlying undirected graph ≈ longest propagation path
    try:
        diameter = nx.diameter(graph.to_undirected())
    except nx.NetworkXError:
        # Disconnected graph — use largest component's diameter
        largest = max(nx.connected_components(graph.to_undirected()), key=len)
        subgraph = graph.subgraph(largest).to_undirected()
        diameter = nx.diameter(subgraph)
    
    # 3× diameter is generous; log factor adds safety for dense cycle interactions
    return max(10, int(3 * diameter * math.log1p(graph.number_of_nodes())))

def run_cmpa(data, graph,tol=1e-6, damping=0.85):
    """
    Runs the CMPA algorithm on a networkx graph, compatible with cyclic graphs.
    
    Uses iterative convergence: repeatedly updates all node scores until they
    stabilize (delta < tol) or max_iter is reached. Cyclic nodes will reach
    a fixed-point equilibrium rather than causing an infinite loop.

    Arguments:
    - data:     dict of {gene_symbol: log_FC score}
    - graph:    networkx MultiDiGraph with 'protein' on nodes, 'type' on edges
    - tol:      convergence tolerance — stops when max score delta < tol (default: 1e-6)
    - damping:  damping factor for cyclic nodes to prevent score explosion (default: 0.85)

    Returns:
    - pandas DataFrame with columns ['hub', 'score']
    """
    is_dag = nx.is_directed_acyclic_graph(graph)
    if not is_dag:
        cycles = list(nx.simple_cycles(graph))
        print(f"Warning: Graph contains {len(cycles)} cycle(s). Using iterative convergence.")
        #print(f"  Cycles detected: {cycles}")
    else:
        print("Graph is a DAG. Running standard topological CMPA.")

    # ── Step 1: Initialise scores ──────────────────────────────────────────────
    for node, attr in graph.nodes(data=True):
        node_protein = attr.get('protein')
        base = data[node_protein] if node_protein in data else 1.0
        graph.nodes[node]['base_score'] = base   # immutable seed value
        graph.nodes[node]['IF_score']   = base   # will be updated each iteration

    # ── Step 2: Determine processing order ────────────────────────────────────
    if is_dag:
        # Topological order guarantees each predecessor is processed before its successor
        ordered_nodes = list(nx.topological_sort(graph))
    else:
        # For cyclic graphs: isolate cycle nodes so we can apply damping to them
        cycle_nodes = set(n for cycle in nx.simple_cycles(graph) for n in cycle)
        # Topologically sort the condensation (SCC DAG), then expand each SCC's members
        # SCC = Strongly Connected Components
        condensation  = nx.condensation(graph)
        ordered_nodes = [
            node
            for scc_idx in nx.topological_sort(condensation)
            for node in condensation.nodes[scc_idx]['members']   # O(1) lookup per SCC
        ]

    # ── Step 3: Iterative convergence ─────────────────────────────────────────
    max_iter = get_max_iter(graph)
    #print('Max iteration:', max_iter)
    for iteration in tqdm(range(max_iter), desc="Iterative convergence"):
        max_delta = 0.0

        for node in ordered_nodes:
            old_score = graph.nodes[node]['IF_score']
            new_score = compute_IF_score(graph, node)

            # Apply damping only to nodes that are part of a cycle
            if not is_dag and node in cycle_nodes:
                new_score = damping * new_score + (1 - damping) * graph.nodes[node]['base_score']

            graph.nodes[node]['IF_score'] = new_score
            max_delta = max(max_delta, abs(new_score - old_score))

        #print(f"Iteration {iteration + 1:>3}: max_delta = {max_delta:.2e}")

        # Early exit once scores have converged
        if max_delta < tol:
            #print(f"Converged after {iteration + 1} iteration(s).")
            break
    #else:
        # print(f"Warning: Reached max_iter={max_iter} without full convergence "
        #       f"(final max_delta={max_delta:.2e}). Scores are approximate.")

    # ── Step 4: Collect results ────────────────────────────────────────────────
    result_list = [
        {'hub': node, 'score': graph.nodes[node]['IF_score']}
        for node in graph.nodes
    ]
    return pd.DataFrame(result_list)

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
    confidence_map = {'Very High':1, 'High': 0.8, 'Medium': 0.6, 'Low':0.3, 'Wrong':0.0}
      
    G = nx.MultiDiGraph()
    
    # add nodes to G
    for node in graph_dict['nodes']:
        #print(node)
        id = node.get('name','')
        label = node.get('properties',{}).get('type',"")
        #label = "".join(word.strip().capitalize() for word in label.split('_'))
        #label=label.capitalize()
        #print(label)
        if label.lower() == "protein":
            protein = node.get('properties',{}).get('name',"")
            protein = protein.split('_')[0]
        else:
            protein = ""
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
            
            subgraph_scores.append(df_subG['score'].sum())
        df_patient_cmpa.iloc[i] = subgraph_scores
    
    # save features file
    #save_path = os.path.join(output_dir, 'ad_subgraph_features.csv')
    df_patient_cmpa.to_csv(output_dir)
    return df_patient_cmpa
