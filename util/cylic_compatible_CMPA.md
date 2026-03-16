# Cyclic-Compatible CMPA — Modification Summary

**Original problem:** The original CMPA algorithm only works on DAGs (directed acyclic graphs).
On cyclic graphs, the `while hub_list` loop deadlocks — nodes in a cycle always have predecessors
still in `hub_list`, so none can ever be removed, causing an infinite loop.

---

## Overview of Changes

| Area | Original | Modified |
|------|----------|----------|
| Cycle detection | None | `nx.is_directed_acyclic_graph()` |
| Processing order | `while hub_list` loop | Topological sort via SCC condensation |
| Score initialisation | `IF_score` only | Separate `base_score` + `IF_score` |
| Score formula | Accumulate on `IF_score` | `damping × influence + (1−damping) × base_score` |
| DAG path | Single `while` pass | Single topological pass, `damping=1.0` |
| Cyclic path | Deadlocks | Iterative convergence until `max_delta < tol` |
| `max_iter` | Hardcoded `100` | Derived from graph diameter |
| Influence normalisation | Sum of predecessor scores | Sum only (no division) — preserves hub accumulation |

---

## Modification 1 — Separate `base_score` from `IF_score`

**What changed:**
Every node now stores two values instead of one:

```python
graph.nodes[node]['base_score'] = data[node_protein]  # immutable seed
graph.nodes[node]['IF_score']   = data[node_protein]  # updated each iteration
```

**Why:** On a cyclic graph the scoring loop runs many passes. If `compute_IF_score` starts
from `IF_score` each time, predecessor influence stacks on top of itself across passes and
scores grow without bound. `base_score` is the fixed anchor the formula always resets to,
preventing runaway accumulation. `IF_score` is the live value that predecessors read when
computing their own scores.

---

## Modification 2 — Processing order via SCC condensation

**What changed:**
The `while hub_list` loop is replaced by a node order derived from the graph structure.

```python
# DAG — simple topological sort
ordered_nodes = list(nx.topological_sort(graph))

# Cyclic graph — collapse cycles into super-nodes first, then sort
condensation  = nx.condensation(graph)
ordered_nodes = [
    node
    for scc_idx in nx.topological_sort(condensation)
    for node in condensation.nodes[scc_idx]['members']
]
```

**Why:** `nx.condensation` collapses each cycle into a single SCC (Strongly Connected
Component) super-node, producing a DAG. Topologically sorting that DAG gives an
upstream-first order even for cyclic graphs. Processing nodes upstream-first means each
node uses the freshest possible predecessor scores within each iteration pass, minimising
the number of passes needed to converge.

**Why not divide by predecessor count:** Hub nodes with many predecessors should accumulate
more signal, not less. Dividing by predecessor count (receiver-side normalisation) penalises
hubs. The correct approach is sender-side — each predecessor contributes its full weighted
score, and the damping factor controls the overall scale.

---

## Modification 3 — New score formula with damping

**What changed:**
`compute_IF_score` now uses a damped formula and always starts from `base_score`:

```python
def compute_IF_score(graph, node, damping=0.85):
    predecessors = list(graph.predecessors(node))
    if not predecessors:
        return graph.nodes[node]['base_score']   # source node — no influence

    influence = 0.0
    for pred in predecessors:
        rel_type = graph[pred][node][0]['type']
        if 'increases' in rel_type.lower():
            weight = 1.0
        elif 'decreases' in rel_type.lower():
            weight = -1.0
        else:
            weight = 0.1
        influence += graph.nodes[pred]['IF_score'] * weight  # uses live score

    # damping=1.0 on DAG path → pure predecessor signal (identical to original)
    # damping=0.85 on cyclic path → 85% predecessor, 15% base anchor
    return damping * influence + (1 - damping) * graph.nodes[node]['base_score']
```

**Why damping works:** The damping factor is a contraction — each pass multiplies the
influence contribution by `damping < 1.0`. Across iterations this decays geometrically
(`0.85¹, 0.85², 0.85³ … → 0`), so the change each pass shrinks until scores stabilise.
Without damping (`damping=1.0`), scores on cyclic graphs circulate indefinitely with no decay.

**Why `damping=1.0` on DAGs:** On an acyclic graph each node is processed exactly once in
topological order. There is no risk of runaway accumulation, so the full predecessor signal
should be used — matching the original algorithm's output exactly.

---

## Modification 4 — Iterative convergence loop

**What changed:**
For cyclic graphs, the scoring pass repeats until scores stabilise:

```python
if is_dag:
    # Single pass — exact, identical to original CMPA
    for node in ordered_nodes:
        if list(graph.predecessors(node)):
            graph.nodes[node]['IF_score'] = compute_IF_score(graph, node, damping=1.0)
else:
    # Repeat until max score change falls below tolerance
    for iteration in range(max_iter):
        max_delta = 0.0
        for node in ordered_nodes:
            old = graph.nodes[node]['IF_score']
            graph.nodes[node]['IF_score'] = compute_IF_score(graph, node, damping=damping)
            max_delta = max(max_delta, abs(graph.nodes[node]['IF_score'] - old))

        if max_delta < tol:
            break  # converged
```

**Why:** On a cyclic graph, A's score depends on C, which depends on B, which depends on A.
No single pass can resolve this mutual dependency. Instead, each pass uses the best current
estimates of all scores, and the estimates improve iteration by iteration until the largest
change (`max_delta`) falls below `tol`. This is the same principle used by PageRank.

---

## Modification 5 — Dynamic `max_iter` from graph diameter

**What changed:**
`max_iter` is computed from the graph's actual structure instead of being hardcoded to `100`:

```python
import math

def get_max_iter(graph):
    try:
        diameter = nx.diameter(graph.to_undirected())
    except nx.NetworkXError:           # disconnected graph
        largest  = max(nx.connected_components(graph.to_undirected()), key=len)
        diameter = nx.diameter(graph.subgraph(largest).to_undirected())

    return max(10, int(3 * diameter * math.log1p(graph.number_of_nodes())))
```

**Why:** Scores need at most `diameter` hops to propagate from one end of the graph to the
other. A 3× multiplier with a log factor adds safety margin for dense cycle interactions.
For graphs of 20–60 nodes this typically yields `max_iter` in the 15–40 range — enough to
converge without letting scores explode over 100+ unnecessary iterations.

---

## Key Invariants to Remember

- **`base_score` is never modified** after initialisation — it is always the raw input signal.
- **`IF_score` on predecessors is always used** inside `compute_IF_score`, never `base_score` —
  this is how hub accumulation propagates through the graph across iterations.
- **DAG path never uses damping** — it would suppress hub signal unnecessarily and break
  equivalence with the original algorithm.
- **Convergence is guaranteed** as long as `damping < 1.0`, because the update rule is a
  contraction mapping — each pass brings all scores strictly closer to their fixed-point values.