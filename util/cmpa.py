# CMPA Algorithm - Cyclic Graph Compatible
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