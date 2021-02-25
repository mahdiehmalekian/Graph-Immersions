"""
Contains some functions for working with dictionaries of format {x: frequency of x}
(function get_minimal is an exception)

@author: Mahdieh Malekian
"""


def add_dict(D, x, k=1):
    """Adds element x to dictionary D k times"""
    if x in D:
        D[x] += k
    else:
        D[x] = k
        
def remove_dict(D, x, k=1):
    """Removes element x from dictionary D k times"""
    D[x] -= k
    
def is_submultiset(L, M):
    """Tells if every entry in L is in M and with >= frequency"""
    for x in L:
        if not x in M or L[x] > M[x]: return False
    return True


def sub_multiset(D, k, already_chosen=None, to_be_ignored=None):
    """"Yields all k-submultisets of a multiset"""
    """Items are chosen from the keys of a dictionary D that are not in 'to_be_ignored'."""
    """D has the form {item: #of its repeatitions}"""
        
    if already_chosen == None: already_chosen = []
    if to_be_ignored == None: to_be_ignored = [] 
    
    if k < 0:
        yield 'none'
        return
    if k == 0:
        yield already_chosen
        return 
    if len(D) == 0:
        yield 'none'
        return   
    for x in D: # Choose an item of D not in to_be_ignored
        if not x in to_be_ignored:
            v = x
            break
    if len(D) - len(to_be_ignored) == 1: # If there is only one item to choose from:
        if D[v] >= k:
            yield already_chosen + [[v, k]]
        else:
            yield 'none'
        return
    
    failed = True # failed will tell if we fail in choosing k items
    # num_use_v is the number of times the item v will be chosen in the choosing of k items
    for num_use_v in [i for i in range(1, min(k, D[v]) + 1)]:   
        for x in sub_multiset(D, k-num_use_v, already_chosen+[[v, num_use_v]], to_be_ignored+[v]):
            if x != 'none':
                failed = False
                yield x
    for x in sub_multiset(D, k, already_chosen, to_be_ignored+[v]):
        if x != 'none':
            failed = False
            yield x
         
    if failed:
        yield 'none'
        
        
        
def pair_up (D, already_paired=None):
    """Yields all the different ways of pairing up items in a dictionary D"""
    """D has the form {item: #of its repetitions}"""
    """Assumes \sum_{d \in D}(#of repetitions of d) is a positive even number"""
    
    if already_paired == None: already_paired = []  
    if D == {}:
        yield already_paired
        return
    for x in D: # Choose an item to be paired
        v = x
        break
        
    if len(D) == 1:
        yield already_paired + [[(v, v), D[v]//2]]
        return
    
    for num_self_paired_v in [i for i in range(D[v]//2 + 1)]: # num_self_paired_v is the number of times v gets paired with itself
        for match in sub_multiset(D, D[v] - 2*num_self_paired_v, already_chosen=[], to_be_ignored =[v]):
            if match != 'none': #update what remains to be paired
                left_over = {}
                visited = [v]
                for x in match:
                    visited.append(x[0])
                    k = D[x[0]] - x[1]
                    if k > 0:
                        left_over [x[0]] = k
                for x in D:
                    if not x in visited:
                        left_over[x] = D[x]
                        
                paired = already_paired + [[(v, x[0]), x[1]] for x in match]
                if num_self_paired_v != 0:
                    paired += [[(v, v), num_self_paired_v]]
                
                for x in pair_up (left_over, paired):
                    yield x
                    
def get_minimal(L):
    """Assumes L is a list of multigraphs with the same underlying simple graph."""
    """Returns the minimal graphs in L"""
    if L == []: return []
    L_num_edge = {}
    minimal, minimal_comb =[], []

    for G in L:
        if G.order() == 0: continue
        m_G = G.num_edge()
        if m_G in L_num_edge:
            if G in L_num_edge[m_G]: continue
            L_num_edge[m_G].append(G)
        else: L_num_edge[m_G] = [G]
    m, M = min(L_num_edge), max(L_num_edge)
    for graph_m in L_num_edge[m]: # every graph in L_num_edge[m] is minimal
        minimal.append(graph_m)
        minimal_comb.append(graph_m.edges())
    
    for num_edg in range(m+1, M+1):
        if not num_edg in L_num_edge: continue
        # L_t is the minimal graphs found in this round, and comb_t is the edge combination of graphs in L_t
        comb_t, L_t = [], []
        # every graph in L_num_edge[num_edg], num_edg > m, is compared against all graphs already in minimal from previous rounds.
        for graph_num in L_num_edge[num_edg]:
            graph_num_is_supergraph = False
            graph_num_comb = graph_num.edges()
            for comb in minimal_comb:
                if is_submultiset(comb, graph_num_comb):
                    graph_num_is_supergraph = True
                    break
            if graph_num_is_supergraph: continue
            L_t.append(graph_num)
            comb_t.append(graph_num_comb)
        minimal += L_t
        minimal_comb += comb_t
    return minimal