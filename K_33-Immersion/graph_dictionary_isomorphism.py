"""
Contains functions for determining if two graphs with dictionary representation, as in
GraphDict or RootedGraphDict classes are isomorphic.

@author: Mahdieh Malekian
"""

import dictionary_operations as dop
 

def unlabeled(G, excluded=None):
    """Returns G[V\excluded] with vertex labels removed"""
    if excluded is None: excluded = []
    L = []
    for v in G:
        if not v in excluded:
            G_v_unlabeled = {}
            for nv in G[v]:
                if not nv in excluded:
                    dop.add_dict(G_v_unlabeled, G[v][nv])
            L.append(G_v_unlabeled)
    return L


def unlabeled_are_equal(G, H, G_excluded=None, H_excluded=None, order_check=None):
    """Tells if by removing vertex labels, G\G_excluded and H\H_excluded are equal"""
    if order_check is None:
        if G_excluded is None: G_excluded = []
        if H_excluded is None: H_excluded = []
        if len(G)-len(G_excluded) != len(H)-len(H_excluded):
            return False
        
    unlabel_G = unlabeled(G, G_excluded)
    unlabel_H = unlabeled(H, H_excluded)
    while unlabel_G:
        unlabel_G_v = unlabel_G.pop()
        if unlabel_G_v not in unlabel_H:
            return False
        unlabel_H.remove(unlabel_G_v)
    return True


def get_most_dominant(G, excluded=None):
    """Returns a vertex with max number of neighbors in G\excluded"""
    if excluded is None:
        excluded = []
    [r, num_good_nbr] = [0, -1]
    for v in G:
        if not v in excluded:
            num_good_nbr_v = 0
            for y in G[v]:
                if not y in excluded:
                    num_good_nbr_v += 1
            if num_good_nbr_v > num_good_nbr:
                r = v
                num_good_nbr = num_good_nbr_v
    return r


def are_locally_similar(G_r, H_s, G_mapped, H_mapped):
    """Tells if 
    1) G_r, H_s restricted to mapped vertices are the same, and
    2) after removing labels of the unmapped vertices, G_r, H_s are the same."""
    G_r_rest = {}
    for nr in G_r:
        if nr in G_mapped:
            nr_index = G_mapped.index(nr)
            if not H_mapped[nr_index] in H_s: return False
            if G_r[nr] != H_s[H_mapped[nr_index]]: return False
        else:
            dop.add_dict(G_r_rest, G_r[nr])
            
    H_s_rest = {}
    for ns in H_s:
        if not ns in H_mapped:
            dop.add_dict(H_s_rest, H_s[ns])
            
    return G_r_rest == H_s_rest

        
def pick_new_roots(G, H, G_roots, H_roots):
    """Yields new roots for G and H which have the same "partition" w.r.t.
    already mapped vertices"""
    G_new_root = get_most_dominant(G, G_roots)
    G_r = G[G_new_root]
    new_root_num_nbr = len(G_r)
    H_candidates = [x for x in H if not x in H_roots and len(H[x]) == new_root_num_nbr]
    if H_candidates == []:
        yield None
        return
    failed = True
    for s in H_candidates:
        H_s = H[s]
        if are_locally_similar(G_r, H_s, G_roots, H_roots):      
            failed = False
            yield (G_new_root, s)
    if failed:
        yield None


def proposed_roots_work(G, H, G_proposed_roots, H_proposed_roots):
    """Tells if
    1) G[G_proposed_roots] is the same as H[H_proposed_roots], and
    2) G, H look similar on the unmapped"""
    if len(G) != len(H):
        return False
    num_roots = len(G_proposed_roots)
    if len(H_proposed_roots) != num_roots:
        return False
    if not unlabeled_are_equal(G, H, order_check=True):
        return False
    for i in range(num_roots):
        if len(G[G_proposed_roots[i]]) != len(H[H_proposed_roots[i]]):
            return False
    for i in range(num_roots):
        if not are_locally_similar(G[G_proposed_roots[i]], H[H_proposed_roots[i]],
                                   G_proposed_roots[:i], H_proposed_roots[:i]):
            return False  
    return unlabeled_are_equal(G, H, G_proposed_roots, H_proposed_roots, order_check=True)

def are_isomorphic(G, H, G_roots=None, H_roots=None, root_check=None, the_same=None):
    """Decides if G, H are isomorphic"""
    if the_same is None and G_roots == H_roots and G == H: return True   
    if G_roots is None: G_roots = []
    if H_roots is None: H_roots = []     
    if not root_check:
        if not proposed_roots_work(G, H, G_roots, H_roots):
            return False          
    if len (G_roots) == len (G): return True   
    for new_roots in pick_new_roots(G, H, G_roots, H_roots):
        if new_roots is None:
            break
        if are_isomorphic(G, H, G_roots+[new_roots[0]], H_roots+[new_roots[1]],
                          root_check=True, the_same=False):
            return True     
    return False


def is_new(G, L):
    """Checks whether G is a "new" graph to L, i.e. no graph isomorphic to G is in L.
     Assumes G has dict type (rather than GraphDict).""" 
    if L == []: return True
    for H in L:
        if G == H or are_isomorphic(G, H):
            return False
    return True


def is_new_r(G, L, root):
    """Checks whether G is a "new" rooted graph to L, i.e. no graph isomorphic to G is in L.
    Assumes G has dict type (rather than RootedGraphDict)."""
    if L == []: return True
    for H in L:
        if G == H or are_isomorphic(G, H, G_roots=[root], H_roots=[root]):
            return False
    return True