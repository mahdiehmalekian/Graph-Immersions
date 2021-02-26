"""
A class that can be used for a dictionary representation of loopless undirected graphs.
Works for simple and multigraphs.

@author: Mahdieh Malekian
"""

import dictionary_operations as dop
import itertools



class GraphDict(dict):
    """A graph_dict is a dictionary representation of graph"""
    """"Assumes the graph is undirected and loopless""" 
    
    def __init__(self, graph_dictionary):
        """Creates a graph dictionary with given graph_dictionary"""
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
        self.graph_dict = graph_dictionary      
        
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
        
    def edge_mult(self, u, v, edge_exists=False):
        """Returns edge multiplicity of (u,v)"""
        return self.graph_dict[v].get(u, 0)
    
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
            self.del_vx(v)
            
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
        
    def e_boundary(self, X):
        """Returns the edge boundary of a set X of vertices in the from {edge:multiplicity}"""
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
        """Adds k edges to the boundary of X in parallel to the existing edges in all possible ways"""
        L = []
        for added in itertools.combinations_with_replacement(self.e_boundary(X), k):
            Gplus = self.copy()
            Gplus.add_edges(added)
            L.append(Gplus)
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
       
    def __str__(self):            
        return str(self.graph_dict)