"""
Contains computations for finding obstructions to rooted W_4 immersion.

@author: Mahdieh Malekian
"""
import datetime
import rooted_graph_dictionary as rgrd
import w4_particulars as w4

def extract_graphs(filepath):
    """Returns the list of rooted graphs stored in M*_rooted.txt files."""
    L = []
    with open(filepath) as file_object:
        next(file_object)
        for line in file_object:
            edges_string, root = line.split(' root:')
            G_root = int(root[0])
            edges_string = edges_string[2:-2]
            G_edges = {}
            pairs_string = edges_string.split("), (")
            for pair in pairs_string:
                edge = tuple(map(int, pair.split(', ')))
                if edge in G_edges:
                    G_edges[edge] += 1
                else:
                    G_edges[edge] = 1
            G = rgrd.make_rooted_graph_dict(edges=G_edges, root=G_root)
            L.append(G)
    return L

M = [] #M[order-5] contains all connected simple rooted graphs of order order
#read the data from text files
for order in [5, 6, 7, 8]:
    filepath = 'connected_simple_rooted_graphs//M'+str(order)+'_rooted.txt'
    M_order = extract_graphs(filepath)
    M.append(M_order)

star_line = '\n\n****************************************************************************************************************************************'

results_file = 'rooted W_4 obstructions.txt'
with open(results_file, 'a') as file_object:
    file_object.write('The following results are used in the proof of Lemma 3.7 of the paper'+
                      '\nDeVos, Malekian, The structure of graphs with no W_4 immersion,'+
                      ' https://arxiv.org/abs/1810.12873.'+star_line+
                      '\nIn the following, all graphs are undirected and loopless and have '+
                      'one distinguished vertex, called the root.'+
                      "\n\nBy a 'well-edge-connected rooted graph', we mean a 3-edge-connected, "+
                      'internally 4-edge-connected rooted multigraph in which deg(root)>=4'+
                      '\nand d(root, root neighbor)>= 5 for any neighbor of the root vertex.'+
                      '\nAlso, if X is a subset of vertices V with |X|>=3 and |V\X|>=3 we have d(X)>=5.'+
                      star_line+
                      '\nComputation for rooted graphs on 5, 6, 7, and 8 vertices:\n'+
                      '(Edge multiplicities are at most eight.)')

time_0 = datetime.datetime.now()

N = [] #N[order-5] contains all graphs in M[order-5] which don't immerse rooted W_4,
O = [] #O[order-5] consists of all "obstructions" arising from the graphs in N[order-5]
exception = [] #exception is any non-Type 2 graph in O
accuracy = 0
for order in [5, 6, 7, 8]: 
    N_order = []
    for G in M[order-5]:
        if w4.has_rw4_im(G, exception=exception, accuracy=accuracy): continue
        N_order.append(G)        
    N.append(N_order)
    O_order = []
    for G in N_order:
        obstruction_G = w4.obstruction(G, exception=exception, accuracy=accuracy)
        if obstruction_G == []: continue
        O_order += obstruction_G        
    for obstruction in O_order:
        if not w4.is_rw4_type2(obstruction, two_econ=True):
            exception.append(obstruction)        
    with open(results_file, 'a') as file_object:
        file_object.write('\n\nThe number of simple connected rooted graphs on '+str(order)+
                          ' vertices without an immersion of rooted W_4 is '+str(len(N_order))+'.'+
                          '\nThe number of well-edge-connected rooted graphs on ' +str(order)+
                          ' vertices without an immersion of rooted W_4 is '+str(len(O_order))+'.')                   
    accuracy = order
    O.append(O_order)
    
time_1 = datetime.datetime.now()
with open(results_file, 'a') as file_object:
    file_object.write('\n\n\nCalculation was done in '+str(time_1-time_0)+'.')


with open(results_file, 'a') as file_object:
    file_object.write(star_line+'\nAfter filtering out Type 2 graphs:\n\n'+
                      '\nIn the family of well-edge-connected rooted graphs, the followings are'+
                      ' the only graphs without an immersion of rooted W_4 which are not Type 2.'+
                      '\n(Note: The edge sets are in the form of {edge: edge_multiplicity})')

for obstruction in exception:
    with open(results_file, 'a') as file_object:
        file_object.write('\n\nThe connected graph with root '
                          +str(obstruction.root)+ ' and edge set '+ str(obstruction.edges()))
        

with open(results_file, 'a') as file_object:
    file_object.write(star_line+
                      '\nThe above graphs are the graphs of Type 3 or Type 4.'+
                      '\nA drawing of these garphs can be found in '+
                      "'drawing rooted W_4 obstructions.ipynb'.")
