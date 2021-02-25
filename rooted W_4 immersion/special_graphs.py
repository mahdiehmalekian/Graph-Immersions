"""
Contains functions which check whether a graph has special types A_k, B_4.
These types are described in 'Forbidding a D_m immersion' section of 
DeVos, Malekian, The structure of graphs with no W4 immersion,
which can be found at https://arxiv.org/pdf/1810.12863.

@author: Mahdieh Malekian
"""
import itertools

def is_type_A_k(G, u, v, k, rest=[], checked_deg_uv=False, checked_deg_rest=False):
    """Takes G and two vertices of it, u,v and decides if there is an arrangement
    x_1x_2 ... x_{|V|-2} of V(G) \{u, v} for which d(u) = d({u, x_1, \ldots, x_i}) =k,
    for every 1<=i<= |V|-2"""
    if not checked_deg_uv:
        if G.deg(u) != k or G.deg(v) != k:
            return [False]
    if not checked_deg_rest:
        for x in G.vertices():
            if x not in [u, v] and G.deg(x) % 2 != 0:
                return [False]
    if rest == []:
        rest = [x for x in G.vertices() if not x in [u, v]]
    m = len(rest)
    not_k_ecut = []
    k_ecut = []
    for perm in itertools.permutations (rest):
        chunk = [u]
        for i in range(m):
            chunk += [perm[i]]
            if chunk in not_k_ecut: break
            if chunk in k_ecut: continue
            if G.e_boundary_size(chunk) == k:
                if i == m-1: return [True, perm]
                else:
                    k_ecut.append(chunk.copy())
            else:
                not_k_ecut.append(chunk.copy())
                break
    return [False]


def is_type_B_4(H, root, root_p, connected=False, two_econ=False):
    """Decides if graph H is type B_4 w.r.t. root, root_p"""
    if not connected:
        if not H.is_connected():return False
    m = H.edge_mult(root, root_p)
    if m == 0:
        mid_6_vx = []
        for root_nbr in H.nbrs(root):
            if H.edge_mult(root, root_nbr) == 3 == H.edge_mult(root_p, root_nbr):
                if H.deg(root_nbr) == 6:
                    mid_6_vx.append(root_nbr)
        if len(mid_6_vx) == 1: H_possible_type = 'B'
        else: return False
    elif m == 1:
        if two_econ or H.is_connected(ignored_edges=[(root, root_p)]):
            H_possible_type = 'A'
        else: return False              
    elif m == 3: H_possible_type = 'C'
    else: return False
    for node in H.vertices():
        node_det_deg = H.detailed_deg(node)
        if node in [root, root_p]:
            if H_possible_type == 'A':
                if node_det_deg != {2:2, 1:1}: return False                 
            elif node_det_deg != {2:1, 3:1}: return False         
        elif node_det_deg != {2:2}:
            if H_possible_type !='B': return False
            if node_det_deg != {3:2}: return False
            if not node in mid_6_vx: return False
    return True
