"""
Contains functions that are used in finding obstructions to immersion of K_3,3.

@author: Mahdieh Malekian
"""

import itertools
import graph_dictionary as grd
import special_graphs as sgr
import dictionary_operations as dop


def is_k33_type1(G, certificate=None):
    """Tells if G has type 1 w.r.t K_3,3 immersion.
    If it is type 1, then with certificate=True we get three other entries, 
    1) first bag of two vertices,
    2) the permutation of other vertices in between,
    3) the last triple of vertices"""
    V = G.vertices()
    n = len(V)
    if n == 6:
        for triple in itertools.combinations(V, 3):
            if G.e_boundary_size(triple) == 4:
                if certificate:
                    return [True, triple, [x for x in V if not x in triple]]
                else:
                    return True
        return False
    even_deg = [node for node in V if G.deg(node)%2 == 0]
    if len(even_deg) < n-6: return False
    for out_of_bags in itertools.combinations(even_deg, n-6):
        in_bags = [v for v in V if not v in out_of_bags]      
        for bag_1 in itertools.combinations(in_bags, 3):
            if G.e_boundary_size(bag_1) != 4: continue
            bag_2 = [v for v in in_bags if not v in bag_1]
            if G.e_boundary_size(bag_2) != 4: continue
            root, root_p = bag_1[0], bag_2[0]
            H = G.copy()
            H.merge(bag_1)
            H.merge(bag_2)
            is_nested_4ecut = sgr.is_type_A_k(H, root, root_p, k=4, rest=out_of_bags,
                                              checked_deg_uv=True, checked_deg_rest=True)
            if is_nested_4ecut[0]:
                if certificate:
                    return [True, bag_1, is_nested_4ecut[1], bag_2]
                else: return True
    return False
       
def is_k33_type2(G, two_econ=None):
    """Tells if G has Type 2 w.r.t. K_3,3 immersion"""
    if two_econ is None:
        if not G.is_connected(): return False    
    for node in G.vertices():
        if G.num_nbrs(node) > 4: return False
    deg_odd, deg_6, deg_4_22, deg_4_else = [], [], [], [] 
    for node in G.vertices():
        node_deg = G.deg(node)
        if node_deg % 2 == 1: deg_odd.append(node)
        elif node_deg == 4:
            detailed_deg_node = G.detailed_deg(node)
            if detailed_deg_node == {2:2}: deg_4_22.append(node)
            else: deg_4_else.append(node)
        elif node_deg == 6: deg_6.append(node)
    if deg_odd == []: return False
    n = G.order()
    num_22 = len(deg_4_22)
    if num_22 < n-8: return False
    if num_22 == n-8: possible_type = ['A']
    elif num_22 == n-7: possible_type = ['A', 'B']
    else: possible_type = ['A', 'B', 'C']
    num_4_else = len(deg_4_else)
    if num_22 + num_4_else < n-5: return False
    if num_22 + num_4_else == n-5:
        if 'A' in possible_type: possible_type.remove('A')
        if 'C' in possible_type: possible_type.remove('C')     
    if deg_6 == [] and 'B' in possible_type: possible_type.remove('B')
    if possible_type == []: return False
    
    bag_cand = [] # Choose candidates for a "bag"
    for odd_node in deg_odd:
        for edge in G.incident_edges(odd_node):
            if G.e_boundary_size(edge) == 5:
                bag_cand.append(edge)
    if bag_cand == []: return False
    
    bags = itertools.combinations(bag_cand, 2)
    for [bag_1, bag_2] in bags:
        bag_1_1, bag_1_2 = bag_1[0], bag_1[1]
        if bag_1_1 in bag_2 or bag_1_2 in bag_2: continue
        bag_2_1, bag_2_2 = bag_2[0], bag_2[1]
        has_chance = False
        if not 'B' in possible_type:
            bag_1_nbrs = G.nbrs(bag_1_1).union(G.nbrs(bag_1_2))
            if bag_2_1 in bag_1_nbrs or bag_2_2 in bag_1_nbrs:
                has_chance = True
        if has_chance or 'B' in possible_type:
            H = G.copy()
            H.merge(bag_1)
            H.merge(bag_2)
            if sgr.is_type_B_4(H, bag_1_1, bag_2_1, connected=True, two_econ=two_econ):
                return True
    return False

def has_enough_high_deg(G):
    """Returns one or two entries. First entry returns False if G after
    k33_reduction will not have enough vertices of appropriate degrees.
    It may return True generously."""
    """If first entry is True, second entry returns the vertices with deg 1"""
    vertices = G.vertices()
    G_deg = {node: G.deg(node) for node in vertices}
    n = len(vertices)
    deg_1 = [node for node in G_deg if G_deg[node] == 1]
    num_deg_1 = len(deg_1)
    num_deg_2 = len({node for node in G_deg if G_deg[node] == 2})
    if n-num_deg_1-num_deg_2 < 6: return [False]
    num_eventually_removed = 0
    for node in vertices:
        if G_deg[node] < 3: continue
        to_be_removed_nbrs = 0
        for node_nbr in G.nbrs(node):
            if G_deg[node_nbr] == 1:
                to_be_removed_nbrs += 1
        if G_deg[node]-to_be_removed_nbrs < 3:
            num_eventually_removed += 1
    if n-num_deg_1-num_deg_2-num_eventually_removed < 6: return [False]
    return [True, deg_1]
    

def k33_reduction(G):
    """Takes a graph G and reduces its 1, 2- edge-cuts, and internal 3-edge-cuts.
    If G after reduction won't have >= 6 vertices, it may simply return an edge."""
    small_graph = grd.GraphDict(graph_dictionary={0:{1:1}, 1:{0:1}})
    high_deg_check = has_enough_high_deg(G)
    if not high_deg_check[0]: return small_graph    
    H = G.copy()
    H.del_vtcs(high_deg_check[1])
    n = H.order()   
    # Reduce 1, 2-edge-cuts
    for j in [1, 2]:
        for k in range(n//2, 0, -1):
            while H.order() >= max([6, 2*k]):
                [min_d_k, side] = H.min_d_k_tuple(k)
                if min_d_k == j:
                    if j == 1: H.one_ecut_rdn_X(side)
                    else: H.two_ecut_rdn_X(side)                                            
                else: break
            if H.order() < 6: return small_graph
        if H.order() < 6: return small_graph
    # Reduce internal 3-edge-cuts
    for k in range(n//2, 1, -1):
        while H.order() >= max([6, 2*k]):
            [min_d_k, side] = H.min_d_k_tuple(k)
            if min_d_k == 3: H.merge(side)
            else: break
        if H.order() < 6: return small_graph
    return H


def is_well_econ(G, reduced=None, m_tuples_checked=None):
    """Checks if graph G
    1) has min degree three, and
    2) is internally 4-edge connected, and
    3) there is no (<=4)-edge-cut with at least three vertices on both sides.
    'm_tuples_checked' is the maximum m for which the edge-boundary of all
    m-tuples of G meet the required size for well-edge-connectivity.
    """
    n = G.order()
    if n < 6: return False
    if reduced and m_tuples_checked is None:
        m_tuples_checked = 2
    if m_tuples_checked is None or m_tuples_checked < 1:
        for node in G.vertices():
            if G.deg(node) < 3: return False
        m_tuples_checked = 1
    if m_tuples_checked < 2:
        for e in G.edge_set():
            if G.e_boundary_size(e) < 4: return False
        m_tuples_checked = 2
    
    min_check = max(3, m_tuples_checked+1)
    max_check = n//2
    if min_check > max_check: return True
    for k in range(min_check, 1+max_check, 1):
        if not G.min_d_k_tuple(k, min_needed=5): return False
    return True    


def has_k33_sbg_six(G, order_check=None):
    """Determines if G on six vertices has a subgraph of K_3,3"""
    if order_check is None and G.order() != 6: return False
    H = G.underlying_simple_complement()
    H_comps_sizes = H.components_sizes()
    for x in H_comps_sizes:
        if x > 3: return False
    if 1 not in H_comps_sizes and 3 not in H_comps_sizes: return False
    return True


def has_k33_im_six(G, reduced=None, well_econ=None, type1=None, type2=None,
                   m_tuples_checked=None):
    """Assumes G has six vertices and decides whether it immerses K_3,3.
    The optional arguments tell if we already know G is reduced,
    or well-edge-connected, or has type 1 or 2.
    'm_tuples_checked' is the maximum m for which the edge-boundary of all
    m-tuples of G meet the required size for well-edge-connectivity.
    """
    if type1 or type2:
        return False
    if has_k33_sbg_six(G): return True
    if well_econ is None:
        if reduced is None:
            if k33_reduction(G).order() != 6: return False
        if not is_well_econ(G, reduced=True, m_tuples_checked=m_tuples_checked):
            return False
    elif not well_econ: return False
    if type2 is None:
        if is_k33_type2(G, two_econ=True): return False

    #high_deg consits of vertices of degree >=5
    high_deg = []
    for node in G.vertices():
        node_deg = G.deg(node)
        if node_deg >= 5: high_deg.append((node, (node_deg-3)//2))
    if high_deg == []: return False
    previous_graphs = [G]
    for (v, num_split) in high_deg:
        v_splits = []
        # v_splits contains graphs obtained by doing at most num_split splits at
        #v at all graphs in previous_graphs
        for j in range(num_split):
            v_j_splits = []
            for g in previous_graphs:
                for g_split in g.split(v, num_split=1):
                    if g_split is None: continue
                    if g_split.order() != 6: continue
                    if g_split in v_j_splits: continue
                    if has_k33_sbg_six(g_split, order_check=True): return True
                    if not is_well_econ(g_split): continue
                    v_j_splits.append(g_split)
            previous_graphs = grd.cut_iso(v_j_splits)
            v_splits += previous_graphs
        previous_graphs = grd.cut_iso(v_splits)
    return False



def has_k33_im(G, reduced=None, well_econ=None, exception=None, accuracy=0,
               m_tuples_checked=None):
    """Determines if G of any order immerses K_3,3 using the information we have.
    As a last resort, it uses has_k33_im_six or has_im_split"""
    """well_econ can take True, False, or None-- meaning that
    G is/is not/we don't know if is well-edge-connected, respectively.
    Similarly, for 'reduced'.
    'exception': the list of known well-edge-connected exceptions,
    'accuracy': the largest n for which all the exceptions of order <=n are known.
    'm_tuples_checked' is the maximum m for which the edge-boundary of all
    m-tuples of G meet the required size for well-edge-connectivity.
    """
    
    try:
        if G.order() < 6: return False
    except:
        return False
    if well_econ or reduced: H = G
    else: H = k33_reduction(G)
    if H == G:
        n = G.order()
        H_m_tuples_checked = m_tuples_checked
    else:
        n = H.order()
        if n < 6: return False
        H_m_tuples_checked = None
      
    if H != G or not well_econ:
        if is_k33_type1(H):
            return False
    if is_k33_type2(H, two_econ=True): return False
    
    if n <= accuracy and exception != None:
        if H == G and well_econ != None:
            H_well_econ = well_econ
        else:
            H_well_econ = is_well_econ(H, reduced=True,
                                       m_tuples_checked=H_m_tuples_checked)
        if H_well_econ:
            for g in exception:
                if H.is_isomorphic(g): return False
            return True
        elif n == 6: return False
        
    if n == 6: return has_k33_im_six(H, reduced=True, type1=False, type2=False,
                                     m_tuples_checked=H_m_tuples_checked)        
    
    return has_im_split(H, exception=exception, accuracy=accuracy)


def has_im_split(G, exception=None, accuracy=None):
    """determines if a graph on at least seven vertices immerses K_3,3 by
    splitting its vertices completely (one at a time), and check if the
    resulting graph on fewer vertices immerses K_3,3"""
    for node in G.vertices():
        for g in G.split(node, G.deg(node)//2):
            if has_k33_im(g, exception=exception, accuracy=accuracy):
                return True
    return False


def is_loser(G, reduced=None, well_econ=None, exception=None, accuracy=0,
             info_underlying=None, m_tuples_checked=None):
    """Determines if G fails in immersing K_3,3"""
    """Uses information about whether other graphs with the same underlying
    simple graph immerse K_3,3.
    'm_tuples_checked' is the maximum m for which the edge-boundary of all
    m-tuples of G meet the required size for well-edge-connectivity.
    """
    n = G.order()
    if n < 6: return True
    if info_underlying is None: info_underlying = [[], [], [], []]
    [wincomb, failcomb, non_iso_winners, non_iso_losers] = info_underlying
    if accuracy >= n and exception != None:
        if well_econ is None:
            well_econ = is_well_econ(G, m_tuples_checked=m_tuples_checked)
            if well_econ: reduced = True
        if well_econ:
            return not has_k33_im(G, reduced=True, well_econ=True,
                                  exception=exception, accuracy=accuracy,
                                  m_tuples_checked=m_tuples_checked)
    G_comb = G.edges()
    for winning_comb in wincomb:
        if dop.is_submultiset(winning_comb, G_comb): return False
    for failing_comb in failcomb:
        if dop.is_submultiset(G_comb, failing_comb): return True
    for winner in non_iso_winners:
        if G.is_isomorphic(winner):
            wincomb.append(G_comb)
            return False
    for loser in non_iso_losers:
        if G.is_isomorphic(loser):
            failcomb.append(G_comb)
            return True
    if has_k33_im(G, reduced=reduced, well_econ=well_econ,
                  exception=exception, accuracy=accuracy,
                  m_tuples_checked=m_tuples_checked):
        non_iso_winners.append(G)
        wincomb.append(G_comb)
        return False
    else:
        non_iso_losers.append(G)
        failcomb.append(G_comb)
        return True


def repair_d_m_tuple(G, m, exception=None, accuracy=0, info_underlying=None,
                     m_tuples_checked=None):
    """Adds edges in parallel to the existing edges in the edge boundary of
    m-tuples in G, so that
    1) if m=1, d(tuple)>=3, and
    2) if m=2, d(tuple)>= 4, and
    3) if m>=3, d(tuple)>=5"""
    """info_underlying is the information we have so far about what combination
    of edge multiplicities make G immerse or not immerse K_3,3.
    'm_tuples_checked' is the maximum m for which the edge-boundary of all
    m-tuples of G meet the required size for well-edge-connectivity.
    """
    if info_underlying is None: info_underlying = [[], [], [], []]    
    L = [G]
    n = G.order()
    if m == 1: needed = 3
    elif m == 2: needed = 4
    elif m >= 3: needed = 5
    for m_tuple in itertools.combinations(G.vertices(), m):
        Lt = []
        for g in L:
            if g in Lt: continue
            dtuple = g.e_boundary_size(m_tuple)
            if dtuple >= needed: Lt.append(g)
            else:
                for g_plus in g.add_dX(m_tuple, needed- dtuple):
                    if g_plus in Lt: continue
                    if n == 6: Lt.append(g_plus)
                    elif is_loser(g_plus, exception=exception, accuracy=accuracy,
                                  info_underlying=info_underlying,
                                  m_tuples_checked=m_tuples_checked):
                        Lt.append(g_plus)
        L = dop.get_minimal(Lt)
    return grd.cut_iso(L)


def repair(G, exception=None, accuracy=0, info_underlying=None):
    """Takes a simple graph G as input, and outputs all edge-minimal graphs
    (up to isomorphism) whose undelying simple graph is G, are 3-edge-connected,
    internally 4-edge-connected, and d(X) >= 5 for X\subset V(G) with
    |X|, |V(G)\X|>= 3, and have edge-multiplicity up to nine, and do not immerse K_3,3"""
    if info_underlying is None: info_underlying = [[], [], [], []]
    L = [G]
    n = G.order()
    for m in range(1, n//2+1):
        Lt = []
        for minimal in L:
            for g in repair_d_m_tuple(minimal, m, exception=exception,
                                      accuracy=accuracy,
                                      info_underlying=info_underlying,
                                      m_tuples_checked=m-1):
                if not g in Lt: Lt.append(g)
        if Lt == []: return []
        Lt = dop.get_minimal(Lt)
        L = grd.cut_iso(Lt)
    if n == 6:
        L = [x for x in L
             if is_loser(x, reduced=True, well_econ=True, exception=exception,
                         accuracy=accuracy, info_underlying=info_underlying)]
    return L

def obstruction(G, exception=None, accuracy=0):
    """Takes a simple graph G as input, and outputs all graphs (up to isomorphism)
    whose undelying simple graph is G, are 3-edge-connected, internally 4-edge-connected,
    and d(X) >= 5 for X\subset V(H) with |X|, |V(H) \X| >= 3,
    and have edge-multiplicity up to nine, and do not immerse K_3,3"""
    H = G.copy()
    info_underlying = [[], [], [], []]
    losers = repair(H, exception=exception, accuracy=accuracy,
                    info_underlying=info_underlying)
    if losers == []: return []
    G_obstruction = []
    for x in grd.cut_iso(losers):
        G_obstruction.append(x)
    for step in range(8*H.num_edge()):
        Lt = []
        for g in losers:
            for e in g.edge_set():
                if g.edge_mult(e[0], e[1]) < 9 :
                    gplus = g.copy()
                    gplus.add_edge(e[0], e[1])
                    if not gplus in Lt: Lt.append(gplus)
        losers_step = [x for x in Lt
                       if is_loser(x, reduced=True, well_econ=True,
                                   exception=exception, accuracy=accuracy,
                                   info_underlying=info_underlying)]
        if losers_step == []: break
        losers = grd.cut_iso(losers_step)
        for x in losers:
            if not x in G_obstruction: G_obstruction.append(x)
    return grd.cut_iso(G_obstruction)