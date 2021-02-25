"""
A class that can be used for a dictionary representation of loopless graphs.
Works for simple and multigraphs.

@author: Mahdieh Malekian
"""

import graph_dictionary_isomorphism as iso
import dictionary_operations as dop
import itertools



class GraphDict(dict):
    """GraphDict uses graph_dictionary for a dictionary representation of graph"""
    """Assumes the graph is undirected and loopless""" 
    """The dictionaries are assumed to be in the form of
    {v: {v_nbr: edge_multiplicity of (v, v_nbr) for v_nbr a neighbor of v}
         for v a vertex of graph}
    
    Example 1:
        The dictionary representation of the graph with vertex set {0, 1, 2, 3, 4}
        and edges [(0, 1), (0, 2), (0, 3)] is
        {0: {1: 1, 2: 1, 3: 1}, 1: {0: 1}, 2: {0: 1}, 3: {0: 1}, 4: {}}.
        
    Example 2:
        The dictionary representation of the graph with vertex set {0, 1, 2, 3}
        and edges [(0, 1), (0, 2), (0, 2), (1, 3), (1, 3), (1, 3), (1, 3), (2, 3), (2, 3)]
        is {0: {1: 1, 2: 2}, 1: {0: 1, 3: 4}, 2: {0: 2, 3: 2}, 3: {1: 4, 2: 2}}."""
    
    def __init__(self, graph_dictionary):
        """Creates a graph dictionary with given graph_dictionary"""
        self.graph_dict = graph_dictionary
        
    def dictionary(self):
        return self.graph_dict       
        
    def vertices(self):
        """returns vertices"""
        return {x for x in self.graph_dict}    
    
    def order(self):
        return len(self.graph_dict)
    
    def edges(self):
        """Returns the edges in the form edge:edge multiplicity"""
        """edges are returned as tuples"""
        """Assumes vertex names are comparable."""
        E = {}
        for x in self.graph_dict:
            for y in self.graph_dict[x]:
                if x < y:
                    E[(x, y)] = self.graph_dict[x][y]
        return E
      
    def num_edge(self):
        """Returns the number of edges."""
        """Assumes vertex names are comparable."""
        count = 0
        for v in self.graph_dict:
            for nv in self.graph_dict[v]:
                if v < nv:
                    count += self.graph_dict[v][nv]
        return count    
    
    def deg(self, v):
        """Returns degree of vertex v of self."""
        count = 0
        for u in self.graph_dict[v]:
            count += self.graph_dict[v][u]
        return count
    
    def detailed_deg(self, v):
        dv_unlabeled = {}
        for u in self.graph_dict[v]:
            dop.add_dict(dv_unlabeled, self.graph_dict[v][u])
        return dv_unlabeled
        

    def edge_set(self):
        """yields the edge set."""
        """"Only one edge per a parallel class"""
        """edges are returned as tuples"""
        """"Assumes vertex names are comparable."""
        for x in self.graph_dict:
            for y in self.graph_dict[x]:
                if x < y:
                    yield (x, y)                    
                    
    def nbrs(self, v):
        """Returns the neighbours of v in self"""
        return {nv for nv in self.graph_dict[v]}
    
    def num_nbrs(self, v):
        """Returns the number of neighbours of v in self"""
        return len(self.graph_dict[v])
        
    def incident_edges(self, v):
        """Returns edges incident with vertex v"""
        return [(v, nv) for nv in self.graph_dict[v]]
              
    def edge_mult(self, u, v, edge_exists=None):
        """Returns edge multiplicity of (u, v)"""
        if edge_exists: return self.graph_dict[v][u]
        return self.graph_dict[v].get(u, 0)
   
    def is_isolated(self, v):
        """Returns True if v is isolated, False otherwise."""
        return len(self.graph_dict[v]) == 0
    
    def underlying_simple(self):
        """Returns the underlying simple graph of self"""
        self_underlying = {}
        for v in self.graph_dict:
            self_underlying[v]={nv:1 for nv in self.graph_dict[v]}
        return GraphDict(self_underlying)
    
    def underlying_simple_complement(self):
        """Returns the complement of the underlying simple graph of self"""
        self_complement = {}
        V_self = self.vertices()
        for v in V_self:
            self_complement[v] = {w:1 for w in V_self 
                           if w!= v and not w in self.graph_dict[v]}
        return GraphDict(self_complement)
    
    def merge(self, L):
        """Returns the loopless graph obtained by merging vertices in L."""
        """The first vertex in L represents all vertices in L in the merged graph"""
        """L is a list of vertces of self."""       
        l, Gl = L[0], {}
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
    
    def add_edge(self, u, v):
        """Adds an edge (u, v) to self"""
        if not u in self.graph_dict:
            self.graph_dict[u] = {}
        if not v in self.graph_dict:
            self.graph_dict[v] = {}
        dop.add_dict(self.graph_dict[u], v)
        dop.add_dict(self.graph_dict[v], u)
            
    def add_edges(self, L):
        """Add all edges in L to self"""
        for e in L:
            self.add_edge(e[0], e[1])            
            
    def del_vx(self, u):
        """Deletes the vertex u from self"""
        for nu in self.graph_dict[u]:
            del self.graph_dict[nu][u]
        del self.graph_dict[u]
        
    def del_vtcs(self, L):
        """Deletes all vertices in L from self"""
        for v in L:
            self.del_vx (v)
            
    def del_edge(self, u, v):
        """Deletes one copy of the edge (u, v) from self"""
        mult_uv = self.graph_dict[u][v]
        if mult_uv == 1:
            del self.graph_dict[u][v]
            del self.graph_dict[v][u]
        else:
            self.graph_dict[u][v] = mult_uv - 1
            self.graph_dict[v][u] = mult_uv - 1
                        
    def del_edges(self, L):
        """Deletes all edges in L from self"""
        for e in L:
            self.del_edge (e[0], e[1])
            
    def is_isomorphic(self, other):
        if self.order() != other.order(): return False
        if self.num_edge() != other.num_edge(): return False
        return iso.are_isomorphic(self.graph_dict, other.graph_dict)
    
    def is_new(self, L):
        """Checks whether self is a "new" graph to L, i.e. no graph isomorphic
        to self is in L"""
        if L == []: return True
        for H in L:
            if self.is_isomorphic(H): return False
        return True

            
    def split(self, v, num_split):
        """Yields all loopless graphs resulting from self by performing
        num_split splits at v"""
        """followed by deleting v in case v ends up with degree <2."""
        """Does not change self."""
        L = []
        D = self.graph_dict[v]
        for comb in dop.sub_multiset(D, 2*num_split):
            for edges_pairing in dop.pair_up ({x[0]: x[1] for x in comb}):
                H_dict = {}
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
                            if H_dict[nv] == {}: del H_dict[nv]
                        del H_dict[v]               
                if iso.is_new(H_dict, L):
                    L.append(H_dict)
                    yield GraphDict(graph_dictionary=H_dict)                    
        if L == []: yield
        
    def e_boundary(self, X):
        """Returns the edge boundary of a set X of vertices in the from
        {edge:multiplicity}"""
        return {(v, u): self.graph_dict[v][u]
                for v in X for u in self.graph_dict[v] if not u in X}

    def e_boundary_size(self, X):
        """Returns the edge boundary size of a set X of vertices"""
        count = 0
        for v in X:
            for u in self.graph_dict[v]:
                if not u in X:
                    count += self.graph_dict[v][u]
        return count
 
    def add_dX(self, X, k):
        """Adds k edges to the boundary of X in parallel to the existing edges
        in all possible ways"""
        L = []
        for added in itertools.combinations_with_replacement(self.e_boundary(X), k):
            Gplus = self.copy()
            Gplus.add_edges(added)
            L.append(Gplus)
            #yield Gplus
        return L
    
    def min_d_k_tuple(self, k, min_needed=None):
        """computes the size of minimum edge-cut separating k-tuples of vertices
         from the rest of the graph"""
        vertices = self.vertices()
        n = len(vertices)    
        k_is_half_order = n == 2*k
        visited = []
        if min_needed == None:
            if k >= n: return
            min_ecut_size = -1
            for k_tuple in itertools.combinations(vertices, k):
                if k_is_half_order and set(k_tuple) in visited:
                    continue
                d_k_tuple = self.e_boundary_size(k_tuple)
                if k_is_half_order:
                    visited.append(vertices.difference(set(k_tuple)))
                if min_ecut_size < 0 or d_k_tuple < min_ecut_size:
                    [min_ecut_size, cert] = [d_k_tuple, k_tuple]
            return [min_ecut_size, cert]
        if k >= n: return False
        for k_tuple in itertools.combinations(vertices, k):
            if k_is_half_order and k_tuple in visited:
                continue
            if self.e_boundary_size(k_tuple) < min_needed:
                return False
            if k_is_half_order:
                visited.append(vertices.difference(set(k_tuple)))
        return True
    
    def dfs(self, node, visited=None, ignored_edges=None):
        """Returns the component of self\ignored_edges containing node (using dfs)"""
        if visited == None: visited = []
        if ignored_edges == None: ignored_edges = []
        if node not in visited:
            visited.append(node)
            for neighbour in self.graph_dict[node]:
                if not (neighbour, node) in ignored_edges:
                    if not (node, neighbour) in ignored_edges:
                        self.dfs(neighbour, visited)
        return visited
    
    def is_connected(self, ignored_edges=None):
        """Tells if self\ignored_edges is connected"""
        for v in self.graph_dict:
            node = v
            break
        return len(self.dfs(node, ignored_edges=ignored_edges)) == self.order()
    
    def components(self):
        """Returns the components of self"""
        comps = []
        remaining = self.vertices()
        while remaining:
            node = remaining.pop()
            node_comp = self.dfs(node)
            for connected_node in node_comp:
                if connected_node in remaining: remaining.remove(connected_node)
            comps.append(node_comp)
        return comps
    
    def components_sizes(self):
        """Returns sizes of the components of self"""
        comps_sizes = {}
        remaining = self.vertices()
        while remaining:
            node = remaining.pop()
            node_comp = self.dfs(node)
            for connected_node in node_comp:
                if connected_node in remaining: remaining.remove(connected_node)
            dop.add_dict(comps_sizes, len(node_comp))
        return comps_sizes
        
    
    def one_ecut_rdn_X(self, X):
        self.del_vtcs(X)
    
    def two_ecut_rdn_X(self, X):
        """Assumes X is a subset of V(self) with d(X)=2."""
        """Returns the graph obtained by replacing X with an edge between the
        endpoints of edge-boundary of X in X^c, if the endpoints are distinct"""
        E_X = self.e_boundary(X)
        if len(E_X) == 2:
            L = []
            for edge in E_X: L.append(edge[1])
            if L[0] != L[1]: self.add_edge(L[0], L[1])
        self.del_vtcs(X)
        
    def __eq__(self, other):
        return self.graph_dict == other.graph_dict
    
    def copy(self):
        self_copy = {}
        for v in self.graph_dict:
            self_copy[v] = self.graph_dict[v].copy()
        return GraphDict(self_copy)
    
    def __hash__(self):
        return hash(repr(self))
       
    def __str__(self):
        return str(self.graph_dict)
    
def make_graph_dict(edges, vertices=None):
    """Takes edges (and vertices) of a graph and returns the GraphDict with the
    given edges (and vertices)"""
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
    return GraphDict(G_dict)
            

def cut_iso(L):
    """"Takes a list L of graphs, and keeps only one of the isomorphic ones"""
    S = []
    for G in L:
        if G.is_new(S):
            S.append(G)
    return S
