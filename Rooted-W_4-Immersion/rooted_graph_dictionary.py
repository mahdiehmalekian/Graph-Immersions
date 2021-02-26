"""
A class that can be used for a dictionary representation of loopless graphs with
one distinguished vertex (called root). It is a child class of GraphDict.
Works for rooted simple and rooted multigraphs.

@author: Mahdieh Malekian
"""

import graph_dictionary as grd
import graph_dictionary_isomorphism as iso
import dictionary_operations as dop

class RootedGraphDict(grd.GraphDict):
    """Deals specifically with rooted graphs"""
    def __init__(self, graph_dictionary, root):
        super().__init__(graph_dictionary)
        self.root = root
        
    def root(self):
        return self.root
    
    def merge(self, L):
        """Returns the loopless graph obtained by merging vertices in L."""
        """The first vertex in L represents all vertices in L in the merged graph
        unless root is in L"""
        """L is a list of vertces of self."""
        if self.root in L: l = self.root
        else: l = L[0]
        Gl = {}
        remaining_vertices = [v for v in self.graph_dict if not v in L]
        
        for v in remaining_vertices:
            Gv = self.graph_dict[v]
            lv_mult = 0
            for x in L:
                if x in Gv:
                    lv_mult += Gv[x]
                    del Gv[x]
            if lv_mult != 0:
                Gv [l] = lv_mult
                Gl [v] = lv_mult
        for v in L: del self.graph_dict[v]
        self.graph_dict[l] = Gl
        
    def is_isomorphic(self, other):
        if self.order() != other.order(): return False
        if self.num_edge() != other.num_edge(): return False
        return iso.are_isomorphic(self.graph_dict, other.graph_dict,
                                  G_roots=[self.root], H_roots=[other.root])
    
    def is_new(self, L):
        if L == []: return True
        for H in L:
            if self == H: return False
        for H in L:
            if self.is_isomorphic(H): return False
        return True
    
    def one_ecut_rdn_X(self, X):
        """Assumes X is a subset of V(self) with d(X)=1."""
        """Returns the graph obtained by removing X (if root not in X)
        or X^c (if root in X)"""
        if not self.root in X:
            self.del_vtcs(X)
        else:
            self.del_vtcs([v for v in self.vertices() if not v in X])
    
    def two_ecut_rdn_X(self, X):
        """Assumes X is a subset of V(self) with d(X)=2."""
        """Returns the graph obtained by replacing X (if root not in X)
        or X^c (if root in X) with an edge between the endpoints of
        edge-boundary of X in X^c, if the endpoints are distinct"""
        if self.root in X:
            side = [v for v in self.vertices() if not v in X]
        else: side = X
        E_side = self.e_boundary(side)
        if len(E_side) == 2:
            L = []
            for edge in E_side: L.append(edge[1])
            if L[0] != L[1]: self.add_edge(L[0], L[1])
        self.del_vtcs(side)
        
    def split(self, v, num_split, min_root_deg=None):
        """Yields all loopless graphs resulting from self by performing num_split
        splits at v, followed by deleting v in case v ends up with degree <2."""
        """Does not change self."""
        L = []
        D = self.graph_dict[v]
        root = self.root
        check_root_needed = False
        if min_root_deg != None:
            deg_root = self.deg(root)
            if v == root:
                if deg_root - 2*num_split < min_root_deg:                    
                    yield
                    return
            else:
                root_v_mult = self.edge_mult(v, root)
                if root_v_mult > 1:
                    max_lost = min(2*num_split, root_v_mult)
                    if deg_root - max_lost < min_root_deg:
                        check_root_needed = True
        for comb in dop.sub_multiset(D, 2*num_split):
            check_root_comb_needed = False
            if check_root_needed:
                max_lost_comb = 0
                for x in comb:
                    if x[0] == root:
                        max_lost_comb = x[1]
                        break
                if max_lost_comb > 1 and deg_root - max_lost_comb < min_root_deg:
                    check_root_comb_needed = True
            for edges_pairing in dop.pair_up({x[0]: x[1] for x in comb}):
                if check_root_comb_needed:
                    final_deg_root = deg_root
                    for element in edges_pairing:
                        [(x, y), xy_mult] = element
                        if x == root and y == root:
                            final_deg_root -= 2*xy_mult
                    if final_deg_root < min_root_deg:
                        continue
                H_dict, root_is_lost = {}, False
                for vertex in self.graph_dict:
                    H_dict[vertex] = self.graph_dict[vertex].copy()
                for element in edges_pairing:
                    [(x, y), xy_mult] = element
                    if x != y: #H.add_edge(x, y) xy_mult times
                        dop.add_dict(H_dict[x], y, k=xy_mult)
                        dop.add_dict(H_dict[y], x, k=xy_mult)                    
                    # H.del_edge(x, v) xy_mult times
                    if H_dict[x][v] == xy_mult: # delete x,v as each other's nbrs
                        del H_dict[x][v]
                        del H_dict[v][x]
                    else:
                        H_dict[x][v] -= xy_mult
                        H_dict[v][x] -= xy_mult
                    #H.del_edge(y, v) k times
                    if H_dict[y][v] == xy_mult: # delete y,v as each other's nbrs
                        del H_dict[y][v]
                        del H_dict[v][y]
                    else:
                        H_dict[y][v] -= xy_mult
                        H_dict[v][y] -= xy_mult              
                if len(H_dict[v]) < 2: #if H.deg(v) < 2, delete v
                    deg_v = 0
                    for nv in H_dict[v]:
                        deg_v += H_dict[v][nv]
                    if deg_v < 2:
                        for nv in H_dict[v]:
                            del H_dict[nv][v]
                            if H_dict[nv] == {}:
                                if nv == root:
                                    root_is_lost = True
                                    break
                                del H_dict[nv]
                        del H_dict[v]            
                if root_is_lost: continue
                if iso.is_new_r(H_dict, L, root):
                    L.append(H_dict)
                    yield RootedGraphDict(root=self.root, graph_dictionary=H_dict)                    
        if L == []: yield

       
    def __eq__(self, other):
        if self.root != other.root: return False
        return self.graph_dict == other.graph_dict
        
    def copy(self):
        self_copy = {}
        for v in self.graph_dict:
            self_copy[v] = self.graph_dict[v].copy()
        return RootedGraphDict(graph_dictionary=self_copy, root=self.root)
    
    def __str__(self):
        return str((self.graph_dict, self.root))
    

def cut_iso(L):
    """"Takes a list L of graphs, and keeps only one of the isomorphic ones"""
    S = []
    for G in L:
        if G.is_new(S):
            S.append(G)
    return S  

def make_rooted_graph_dict(edges, root, vertices=None):
    """Takes edges (and vertices) of a graph and returns the RootedGraphDict
    with the given edges (and vertices)"""
    """Assumes 'edges' is either a dictionary in form {edge:edge_multiplicity} 
    or is a list [(edge)] or a tuple ((edge1)).
    
    Example: The edges [(0, 2), (0, 2), (1, 3)] can be entered just like that, 
    or as a tuple, or as {(0, 2): 2, (1, 3): 1}.
    """
    if vertices is None:
        G_dict = {}
    else:
        G_dict = {node:{} for node in vertices}
    for (x, y) in edges:
        if vertices is None:
            if not x in G_dict: G_dict[x] = {}
            if not y in G_dict: G_dict[y] = {}
        try:
            e_mult = edges[(x, y)]
        except:
            e_mult = 1
        dop.add_dict(G_dict[x], y, k=e_mult)
        dop.add_dict(G_dict[y], x, k=e_mult)
    return RootedGraphDict(G_dict, root=root)