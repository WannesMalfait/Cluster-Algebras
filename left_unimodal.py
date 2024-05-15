# %%

# Compute the elements of \Xi_N
n = 6
perms = [[] for _ in range(2**(n-1))]
for i in range(2**(n-1)):
    perm = [0 for _ in range(n)]
    # Element of Xi_N can be represented as binary number:
    # - 0 indicates go up.
    # - 1 indicates go down.
    bits = i.bit_count()
    min_val = bits + 1
    max_val = min_val
    perm[0] = min_val
    for j in range(1,n):
        if (i >> (j - 1)) & 1 == 0:
            max_val += 1
            perm[j] = max_val
        else:
            min_val -=1
            perm[j] = min_val
    perms[i] = perm

# %%

def adjacent(perm1: list[int], perm2: list[int]) -> bool:
    # Check if the two permutation differ by a permutation.
    for i in range(n-1):
        if perm1[i] == perm2[i]:
            continue
        swapped = perm1[i+1] == perm2[i] and perm1[i] == perm2[i+1]
        if not swapped:
            return False
        # check that all the rest is equal
        for j in range(i+2, n-1):
            if perm1[j] != perm2[j]:
                return False
        return True
    return False      

def is_mutation(perm1: list[int], perm2: list[int], eta: list[int] = [0,1,0,0,1,0]) -> bool:
    for i in range(n-1):
        if perm1[i] == perm2[i]:
            continue
        # Mutation only happens if in the same level set.
        return eta[perm1[i] -1] == eta[perm2[i] -1]
    return False


# %%

import networkx as nx

G = nx.DiGraph()
G.add_nodes_from(range(len(perms)))
for i, perm in enumerate(perms):
    for j, other in enumerate(perms):
        if i >= j:
            continue
        if adjacent(perm, other):
            # Sort lexicographically
            if str(perm) > str(other):
                G.add_edge(i, j)
            else:
                G.add_edge(j, i)
pos = nx.nx_pydot.graphviz_layout(G, prog="dot", root=0)
node_labels = { i : "$"+str(perm)+"$" for i,perm in enumerate(perms)}
edge_styles =  {e : "" if is_mutation(perms[e[0]], perms[e[1]]) else "dashed" for e in G.edges}
# labels = {}
nx.draw_networkx(G,pos=pos, labels=node_labels)



# %%
nx.to_latex(G,pos=pos,node_label=node_labels, edge_options=edge_styles)

# %%
