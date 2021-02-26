"""
Contains computations for finding obstructions to K_3,3 immersion.

@author: Mahdieh Malekian
"""
import datetime
import graph_dictionary as grd
import k33_particulars as k33



def extract_graphs(filepath):
    """Returns the list of graphs stored in M*.txt files."""
    L = []
    with open(filepath) as file_object:
        next(file_object)
        for line in file_object:
            edges_string = line[2:-3]
            G_edges = {}
            pairs_string = edges_string.split("), (")
            for pair in pairs_string:
                edge = tuple(map(int, pair.split(', ')))
                if edge in G_edges:
                    G_edges[edge] += 1
                else:
                    G_edges[edge] = 1
            G = grd.make_graph_dict(edges=G_edges)
            L.append(G)
    return L

M = [] #M[order-6] contains all simple connected graphs of the given order
for order in [6, 7, 8, 9]:
    filepath = 'connected_simple_graphs//M'+str(order)+'.txt'
    M_order = extract_graphs(filepath)
    M.append(M_order)

star_line = '\n\n\n****************************************************************************************************************************************'

results_file = 'K_3,3 obstructions.txt'
with open(results_file, 'a') as file_object:
    file_object.write('The following results are used in the proof of Lemma 3.9 of the paper'+
                      '\nDeVos, Malekian, The structure of graphs with no K_3,3 immersion,'+
                      ' https://arxiv.org/abs/1810.12873.'+star_line+
                      '\nIn the following, all graphs are undirected and loopless.'+
                      "\n\nBy a 'well-edge-connected graph', we mean a 3-edge-connected, "+
                      'internally 4-edge-connected multigraph where for any '+
                      '\nX a subset of vertices V with |X|>=3 and |V\X|>=3 we have d(X)>=5.'+
                      star_line+
                      '\nComputation for graphs on 6, 7, 8, and 9 vertices:\n'+
                      '(Edge multiplicities are at most nine.)')

time_0 = datetime.datetime.now()

N = [] #N[order-6] contains all graphs in M[i-6] which don't immerse K_3,3,
O = [] #O[order-6] consists of all "obstructions" arising from the graphs in N[order-6]
exception = [] #exception is any non-Type 2 graph in O
accuracy = 0
for order in [6, 7, 8, 9]: 
    N_order = []
    for G in M[order-6]:
        if k33.has_k33_im(G, exception=exception, accuracy=accuracy): continue
        N_order.append(G)        
    N.append(N_order)
    O_order = []
    for G in N_order:
        obstruction_G = k33.obstruction(G, exception=exception, accuracy=accuracy)
        if obstruction_G == []: continue
        O_order += obstruction_G        
    for obstruction in O_order:
        if not k33.is_k33_type2(obstruction, two_econ=True):
            exception.append(obstruction)        
    with open(results_file, 'a') as file_object:
        file_object.write('\n\nThe number of simple connected graphs on '+str(order)+
                          ' vertices without an immersion of K_3,3 is '+str(len(N_order))+'.'+
                          '\nThe number of well-edge-connected graphs on ' +str(order)+
                          ' vertices without an immersion of K_3,3 is '+str(len(O_order))+'.')                   
    accuracy = order
    O.append(O_order)
    
time_1 = datetime.datetime.now()
with open(results_file, 'a') as file_object:
    file_object.write('\n\n\nCalculation was done in '+str(time_1-time_0)+'.')


with open(results_file, 'a') as file_object:
    file_object.write(star_line+
                      '\nAfter filtering out Type 2 graphs:\n\n'+
                      '\nIn the family of well-edge-connected graphs, the followings are'+
                      ' the only graphs without an immersion of K_3,3 which are not Type 2.'+
                      '\n(Note: The edge sets are in the form of {edge: edge_multiplicity})')

for obstruction in exception:
    with open(results_file, 'a') as file_object:
        file_object.write('\n\nThe connected graph with edge set '+ str(obstruction.edges()))
        
with open(results_file, 'a') as file_object:
    file_object.write(star_line+
                      '\nThe above graphs are the graphs of Type 3 or Type 4.'+
                      '\nA drawing of these garphs can be found in '+
                      "'drawing K_3,3 obstructions.ipynb'.")
