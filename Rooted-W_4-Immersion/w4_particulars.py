"""
Contains functions that are used in finding obstructions to immersion of rooted W_4.

@author: Mahdieh Malekian
"""

import itertools
import rooted_graph_dictionary as rgrd
import special_graphs as sgr
import dictionary_operations as dop

def is_rw4_type1(G, certificate=None):
    """Determines if G has rooted W_4 Type 1.
    If it is type 1, with certificate=True, we get three other entries:
    1) first bag of two vertices,
    2) the permutation of other vertices in between,
    3) the last bag of three vertices"""
    V = G.vertices()
    n = len(V)
    root = G.root
    if n == 5:
        for edge in G.incident_edges(root):
            if G.e_boundary_size(edge) == 4:
                if certificate:
                    return [True, edge, [], [x for x in V if not x in edge]]
                else:
                    return True
        return False
    even_deg = [node for node in G.vertices() if node!= root and G.deg(node)%2 == 0]
    if len(even_deg) < n-5: return False
    for out_of_bags in itertools.combinations(even_deg, n-5):
        nonroot_in_bags = [v for v in G.vertices() if v!= root and not v in out_of_bags]
        for root_nbr in nonroot_in_bags:
            if G.e_boundary_size([root, root_nbr]) != 4: continue            
            triple = [v for v in nonroot_in_bags if v != root_nbr]
            if G.e_boundary_size(triple) != 4: continue
            H = G.copy()
            H.merge([root, root_nbr])
            H.merge(triple)
            is_nested_4ecut = sgr.is_type_A_k(H, root, triple[0], k=4, rest=out_of_bags,
                                              checked_deg_uv=True, checked_deg_rest=True)
            if is_nested_4ecut[0]:
                if certificate:
                    return [True, [root, root_nbr], is_nested_4ecut[1], triple]
                else:
                    return True
    return False
       
def is_rw4_type2(G, two_econ=None):
    """Tells if G has Type 2 w.r.t. rooted W_4 immersion"""
    root = G.root
    root_det_deg = G.detailed_deg(root)
    if root_det_deg == {2:2, 1:1}: possible_type = ['A', 'C']
    elif root_det_deg == {2:1, 3:1}: possible_type = ['B', 'C']
    else: return False
    
    if two_econ is None:
        if not G.is_connected(): return False
            
    G_detailed_deg = {node: G.detailed_deg(node) for node in G.vertices()}
    n = G.order()
    num_22 = len([node for node in G_detailed_deg if G_detailed_deg[node] == {2:2}])
    if num_22 < n-5: return False
    if num_22 == n-5: 
        if 'C' in possible_type: possible_type.remove('C')
    if possible_type == []: return False
    
    num_211 = len([node for node in G_detailed_deg if G_detailed_deg[node] == {2:1, 1:2}])
    if num_22 + num_211 < n-4: return False
    if num_22 + num_211 == n-4:
        if 'A' in possible_type: possible_type.remove('A')
        if 'C' in possible_type: possible_type.remove('C')
    if possible_type == []: return False
    
    num_deg6 = len([node for node in G_detailed_deg
                    if G_detailed_deg[node] == {3:2}
                    or G_detailed_deg[node] == {3:1, 2:1, 1:1}])
    if num_deg6 == 0 and 'B' in possible_type: possible_type.remove('B')
    if possible_type == []: return False
    
    if G_detailed_deg[root] == {2:2, 1:1}: # Choose vertices that should be in a bag if G is rw4-type 2
        for root_nbr in G.nbrs(root):
            if G.edge_mult(root, root_nbr) == 1:
                in_bag_cand = [root_nbr]
                break
    elif G_detailed_deg[root] == {2:1, 3:1}:
        for root_nbr in G.nbrs(root):
            if G.edge_mult(root, root_nbr) == 3:
                in_bag_cand = [root_nbr_nbr for root_nbr_nbr in G.nbrs(root_nbr)
                if root_nbr_nbr!=root]
                break
    # Choose candidates for a "bag"
    bag_cand = [[in_bag_1, in_bag_2] for in_bag_1 in in_bag_cand
                for in_bag_2 in G.nbrs(in_bag_1)
                if in_bag_2!=root and G.e_boundary_size([in_bag_1, in_bag_2])==5]
    visited = []
    for bag in bag_cand:
        if bag in visited: continue
        H = G.copy()
        H.merge(bag)
        root_p = bag[0]
        if sgr.is_type_B_4(H, root, root_p, connected=True, two_econ=two_econ):
            return True
        visited.append(bag)
        visited.append([bag[1], bag[0]])
    return False


def has_enough_high_deg(G):
    """Returns one or two entries. First entry returns False if G after
    rw_4_reduction will not have enough vertices of appropriate degrees.
    It may return True generously."""
    """If first entry is True, second entry returns the vertices with deg 1."""
    root = G.root
    if G.deg(root) < 4: return [False]
    vertices = G.vertices()
    G_deg = {node: G.deg(node) for node in vertices}
    n = len(vertices)
    deg_1 = [node for node in G_deg if G_deg[node] == 1]
    num_deg_1 = len(deg_1)
    num_deg_2 = len({node for node in G_deg if G_deg[node] == 2})
    if n-num_deg_1-num_deg_2 < 5: return [False]
    num_eventually_removed = 0
    for node in vertices:
        if G_deg[node] < 3: continue
        to_be_removed_nbrs = 0
        for node_nbr in G.nbrs(node):
            if G_deg[node_nbr] == 1:
                to_be_removed_nbrs += 1
        if node == root:
            if G_deg[node]-to_be_removed_nbrs < 4: return [False]
        elif G_deg[node]-to_be_removed_nbrs < 3:
            num_eventually_removed += 1
    if n-num_deg_1-num_deg_2-num_eventually_removed < 5: return [False]
    return [True, deg_1]
    

def rw4_reduction(G):
    """Takes a graph G with root r and if d(r) >= 4, reduces its 1, 2- edge-cuts,
    and internal 3-edge-cuts.
    If G after reduction won't have >= 5 vertices, it may simply return an edge."""
    r = G.root
    small_graph = rgrd.RootedGraphDict(root=r, graph_dictionary={r:{r+1:1}, r+1:{r:1}})
    high_deg_check = has_enough_high_deg(G)
    if not high_deg_check[0]: return small_graph
    
    H = G.copy()
    H.del_vtcs(high_deg_check[1])
    n = H.order()
    # Reduce 1, 2-edge-cuts
    for j in [1, 2]:
        for k in range(n//2, 0, -1):# look at the edge_boundary of k-tuples of vertices
            while H.order() >= max([5, 2*k]):
                [min_d_k, side] = H.min_d_k_tuple(k)
                if min_d_k == j:
                    if j == 1: H.one_ecut_rdn_X(side)
                    else: H.two_ecut_rdn_X(side)                                            
                else: break
            if H.order() < 5: return small_graph
        if H.order() < 5 or H.deg(r) < 4: return small_graph
      
    # Reduce internal 3-edge-cuts
    for k in range(n//2, 1, -1):
        while H.order() >= max([5, 2*k]):
            [min_d_k, side] = H.min_d_k_tuple(k)
            if min_d_k == 3:
                if r in side: side = [v for v in H.vertices() if not v in side]
                H.merge(side)
            else: break
        if H.order() < 5: return small_graph
    return H

def is_well_econ(G, reduced=None, root_checked=None, m_tuples_checked=None):
    """Takes a rooted graph and checks if
    1) d(root)>= 4,
    2) every pair of root and one of its neighbors has edge boundary of size at least five,
    3) graph has min degree three,
    4) graph is internally 4-edge connected,
    5) there is no (<=4)-edge-cut separating at least three vertices opposite others"""
    """If root_checked is True, (1), (2) are skipped.
    'm_tuples_checked' is the maximum m for which the edge-boundary of all
    m-tuples of G meet the required size for well-edge-connectivity.
    """
    n = G.order()
    if n < 5: return False
    root = G.root
    if not root_checked:
        if G.deg(root) < 4: return False
        for root_nbr in G.nbrs(root):
            if G.e_boundary_size([root, root_nbr]) < 5: return False
    if reduced and m_tuples_checked is None:
        m_tuples_checked = 2
    if m_tuples_checked is None or m_tuples_checked < 1:
        for node in G.vertices():
            if node != root and G.deg(node) < 3: return False
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


def has_rw4_sbg_five(G, order_check=None):
    """Determines if G on five vertices has a subgraph of rooted W_4"""
    if not order_check and G.order() != 5: return False
    root = G.root
    if G.num_nbrs(root) != 4: return False
    for node in G.vertices():
        if node != root and G.num_nbrs(node) < 3 : return False
    return True

def has_rw4_im_five(G, reduced=None, well_econ=None, type1=None, type2=None,
                    root_checked=None, m_tuples_checked=None):
    """Assumes G with root root has five vertices and decides whether it has a
    rooted immersion of W_4.
    The optional arguments tell if we already know G is reduced,
    or well-edge-connected, or has Type 1 or Type 2.
    If root_checked is True, it means d(root)>=4 and d(root, root_nbr)>=5, for
    any root_nbr a neighbor of the root vertex.
    'm_tuples_checked' is the maximum m for which the edge-boundary of all
    m-tuples of G meet the required size for well-edge-connectivity. 
    """
    if type1 or type2:
        return False
    elif is_rw4_type1(G):
        return False
    if well_econ is None:
        if reduced is None:
            if rw4_reduction(G).order() != 5: return False
        if not is_well_econ(G, reduced=True, root_checked=root_checked,
                            m_tuples_checked=m_tuples_checked):
            return False
    elif not well_econ: return False
    if has_rw4_sbg_five(G): return True
    if type2 is None:
        if is_rw4_type2(G, two_econ=True): return False
    
    #high_deg consits of all nonroot vertices of degree >=5, and the root vertex if its degree >=6
    high_deg = []
    for node in G.vertices():
        node_deg = G.deg(node)
        if node != G.root:
            if node_deg >= 5: high_deg.append((node, (node_deg-3)//2))
        elif node_deg >= 6: high_deg.append((node, (node_deg-4)//2))
    if high_deg == []: return False
    previous_graphs = [G]
    for (v, num_split) in high_deg:
        v_splits = []
        # v_splits is a list consisting of graphs obtained by doing at most num_split splits at v at all graphs in previous_graphs
        for j in range(num_split):
            v_j_splits = []
            for g in previous_graphs:
                for g_split in g.split(v, num_split=1, min_root_deg=4):
                    if g_split is None: continue
                    if g_split.order() != 5: continue
                    if g_split in v_j_splits: continue
                    if has_rw4_sbg_five(g_split, order_check=True): return True
                    if not is_well_econ(g_split): continue   
                    v_j_splits.append(g_split)
            previous_graphs = rgrd.cut_iso(v_j_splits)
            v_splits += previous_graphs
        previous_graphs = rgrd.cut_iso(v_splits)
    return False


def has_rw4_im(G, reduced=None, well_econ=None, exception=None, accuracy=0,
               root_checked=None, m_tuples_checked=None):
    """Determines if a graph immerses rooted W_4, using the information we have so far.
    As a last resort, it uses has_rw4_im_five or has_im_split"""
    """exception gives the list of the known well-edge-connected exceptions,
    accuracy indicates the largest n for which all the exceptions of order <=n are known."""
    """well_econ can be either True, False, or None-- meaning that G is/is not/
    we don't know if is well-edge-connected, respectively.
    If root_checked is True, it means d(root)>=4 and d(root, root_nbr)>=5, for
    any root_nbr a neighbor of the root vertex.
    'm_tuples_checked' is the maximum m for which the edge-boundary of all
    m-tuples of G meet the required size for well-edge-connectivity.
    """
    try:
        if G.order() < 5: return False
    except:
        return False
    if well_econ or reduced: H = G
    else: H = rw4_reduction(G)      
    if H == G:
        n = G.order()
        H_root_checked, H_m_tuples_checked = root_checked, m_tuples_checked
    else:
        n = H.order()
        if n < 5: return False
        H_root_checked, H_m_tuples_checked = None, None
        
    if H != G or not well_econ:
        if is_rw4_type1(H):
            return False
    if is_rw4_type2(H, two_econ=True): return False
    
    if accuracy >= n and exception != None:
        if H == G and well_econ != None:
            H_well_econ = well_econ
        else:
            H_well_econ = is_well_econ(H, reduced=True,
                                       root_checked=H_root_checked,
                                       m_tuples_checked=H_m_tuples_checked)
        if H_well_econ:
            for g in exception:
                if H.is_isomorphic(g): return False
            return True
        elif n == 5:
            return False
    if n == 5: return has_rw4_im_five(H, reduced=True, type1=False, type2=False,
                                      root_checked=H_root_checked,
                                      m_tuples_checked=H_m_tuples_checked)
       
    return has_im_split(H, exception=exception, accuracy=accuracy)


def has_im_split(G, exception, accuracy):
    """determines if a graph on at least six vertices immerses rooted W_4 by
    splitting its nonroot vertices completely (one at a time), and check if the
    resulting graph on fewer vertices immerses rooted W_4"""
    for node in G.vertices():
        if node == G.root: continue
        for g in G.split(node, num_split=G.deg(node)//2, min_root_deg=4):
            if has_rw4_im(g, exception=exception, accuracy=accuracy): return True
    return False


def is_loser(G, reduced=None, well_econ=None, exception=None, accuracy=0,
             info_underlying=None, root_checked=None, m_tuples_checked=None):
    """Determines if G fails in immersing rooted W_4"""
    """Uses information about whether other graphs with the same underlying
    simple graph immerse W_4.
    If root_checked is True, it means d(root)>=4 and d(root, root_nbr)>=5, for
    any root_nbr a neighbor of the root vertex.
    'm_tuples_checked' is the maximum m for which the edge-boundary of all
    m-tuples of G meet the required size for well-edge-connectivity.
    """
    n = G.order()
    if n < 5: return True
    if info_underlying is None: info_underlying=[[], [], [], []]
    [wincomb, failcomb, non_iso_winners, non_iso_losers] = info_underlying
    if accuracy >= n and exception != None:
        if well_econ is None:
            well_econ = is_well_econ(G, root_checked=root_checked,
                                     m_tuples_checked=m_tuples_checked)
            if well_econ: reduced = True
        if well_econ:
            return not has_rw4_im(G, reduced=True, well_econ=True,
                                  exception=exception, accuracy=accuracy,
                                  root_checked=root_checked,
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
    if has_rw4_im(G, reduced=reduced, well_econ=well_econ,
                  exception=exception, accuracy=accuracy,
                  root_checked=root_checked, m_tuples_checked=m_tuples_checked):
        non_iso_winners.append(G)
        wincomb.append(G_comb)
        return False
    else:
        non_iso_losers.append(G)
        failcomb.append(G_comb)
        return True
    
def repair_root(G, exception=None, accuracy=0, info_underlying=None):
    """Takes a rooted graph G and outputs all edge-minimal graphs
    (up to isomorphism) whose undelying simple graph is the same as G's and
    1) d(root)>= 4, and
    2) every pair of root and one of its neighbors has edge boundary of size
        at least five."""
    r, dr, n = G.root, G.deg(G.root), G.order()
    if info_underlying is None: info_underlying = [[],[],[],[]]
    dr_4 = []
    if dr >= 4: dr_4.append(G)
    else:
        if n == 5:
            dr_4 = [G_plus for G_plus in G.add_dX([r], 4-dr)]         
        else:
            dr_4 = [G_plus for G_plus in G.add_dX([r], 4-dr)
            if is_loser(G_plus, exception=exception, accuracy=accuracy,
                        info_underlying=info_underlying)]
    partially_repaired = rgrd.cut_iso(dr_4)
    for nr in G.nbrs(r):
        d_rnr_5 = []
        for graph in partially_repaired:
            d_rnr = graph.e_boundary_size([r, nr])
            if d_rnr >= 5 and not graph in d_rnr_5:
                d_rnr_5.append(graph)
            else:
                for graph_plus in graph.add_dX([r, nr], 5-d_rnr):
                    if not graph_plus in d_rnr_5:
                        d_rnr_5.append(graph_plus)
        partially_repaired = dop.get_minimal(d_rnr_5)
    L = rgrd.cut_iso(partially_repaired)
    if n != 5:
        L = [x for x in L
             if is_loser(x, exception=exception, accuracy=accuracy,
                         info_underlying=info_underlying, root_checked=True)]
    return L


def repair_d_m_tuple(G, m, exception=None, accuracy=0, info_underlying=None,
                     root_checked=None, m_tuples_checked=None):
    """Adds edges in parallel to the existing edges in the edge boundary of m-tuples in G,
    so that if m=1, d(tuple)>=3, if m=2, d(tuple)>= 4, and if m>=3, d(tuple)>=5"""
    """info_underlying is the information we have so far about what combination
    of edges make G immerse or not immerse rooted W_4.
    If root_checked is True, it means d(root)>=4 and d(root, root_nbr)>=5, for
    any root_nbr a neighbor of the root vertex.
    'm_tuples_checked' is the maximum m for which the edge-boundary of all
    m-tuples of G meet the required size for well-edge-connectivity.
    """
    if info_underlying is None: info_underlying = [[], [], [], []]
    L = [G]
    n = G.order()
    if m == 1:
        needed = 3
        L_m = [[v] for v in G.vertices() if v!=G.root]
    elif m == 2:
        needed = 4
        L_m = [e for e in G.edge_set() if not G.root in e]
    elif m >= 3:
        needed = 5
        L_m = [m_tuple for m_tuple in itertools.combinations(G.vertices(), m)]
    for m_tuple in L_m:
        Lt = []
        for g in L:
            dtuple = g.e_boundary_size(m_tuple)
            if dtuple >= needed: Lt.append(g)
            else:
                for g_plus in g.add_dX(m_tuple, needed-dtuple):
                    if g_plus in Lt: continue
                    if n == 5: Lt.append(g_plus)
                    elif is_loser(g_plus,
                                  exception=exception, accuracy=accuracy,
                                  info_underlying=info_underlying,
                                  root_checked=root_checked,
                                  m_tuples_checked=m_tuples_checked):
                        Lt.append(g_plus)
        L = dop.get_minimal(Lt)
    return rgrd.cut_iso(L)


def repair(G, exception=None, accuracy=0, info_underlying=None):
    """Takes a rooted simple graph G as input, and outputs all edge-minimal graphs
    (up to isomorphism) whose undelying simple graph is G, are 3-edge-connected,
    internally 4-edge-connected, and d(X) >= 5 for X\subset V(G) with
    |X|, |V(G) \X| >= 3, and have edge-multiplicity up to eight,
    and do not immerse rooted W_4"""
    if info_underlying is None: info_underlying = [[], [], [], []]
    L = repair_root(G, exception=exception, accuracy=accuracy,
                    info_underlying=info_underlying)
    n = G.order()
    for m in range(1, n//2 +1):
        Lt = []
        for minimal in L:
            for g in repair_d_m_tuple(minimal, m,
                                      exception=exception, accuracy=accuracy,
                                      info_underlying=info_underlying,
                                      root_checked=True, m_tuples_checked=m-1):
                if g not in Lt: Lt.append(g)
        if Lt == []: return []
        Lt = dop.get_minimal(Lt)
        L = rgrd.cut_iso(Lt)
    if n == 5: L = [x for x in L
                     if is_loser(x, reduced=True, well_econ=True,
                                 exception=exception, accuracy=accuracy,
                                 info_underlying=info_underlying)]
    return L

def obstruction(G, exception=None, accuracy=0):
    """Takes a rooted simple graph G as input, and outputs all graphs (up to isomorphism)
    whose undelying simple graph is G, are 3-edge-connected, internally 4-edge-connected,
    and d(X) >= 5 for X\subset V(H) with |X|, |V(H) \X| >= 3,
    and have edge-multiplicity up to eight, and do not immerse rooted W_4"""
    H = G.copy()
    info_underlying = [[], [], [], []]

    losers = repair(H, exception=exception, accuracy=accuracy,
                    info_underlying=info_underlying)
    if losers == []: return []
    G_obstruction = losers
    for step in range(7*H.num_edge()):
        Lt = []
        for g in losers:
            for e in g.edge_set():
                if g.edge_mult(e[0], e[1]) < 8 :
                    gplus = g.copy()
                    gplus.add_edge(e[0], e[1])
                    Lt.append(gplus)
        losers_step = [x for x in Lt
                       if is_loser(x, reduced=True, well_econ=True,
                                   exception=exception, accuracy=accuracy,
                                   info_underlying=info_underlying)]
        if losers_step == []: break
        losers = rgrd.cut_iso(losers_step)
        for x in losers:
            if not x in G_obstruction: G_obstruction.append(x)
    return rgrd.cut_iso(G_obstruction)