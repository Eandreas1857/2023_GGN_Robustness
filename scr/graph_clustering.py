import matplotlib.pyplot as plt
import networkx as nx
from IPython.display import Image, display
import random

import sys
sys.path.insert(0,'/home/elizabeth/Desktop/GIT/Optimized_Ncut_Directed_and_Undirected/src')
from Clustering_by_weighted_cuts_in_directed_graphs import *


def find_best_clustering(G, start_set, stop_set, network_filename, top_k, nodelist = None, data = None, in_out_degree = 'out', save_file = True):

    if nodelist == None:
        nodelist = list(G.nodes())

    W = asym_weight_matrix(G, nodelist, data)
    D = asym_weighted_degree_matrix(G, nodelist, data, in_out_degree)
    L = Hermitian_normalized_Laplacian(D, W, D)
    Y = k_smallest_eigvec(nodelist, L, 2)

    eigv = sorted(Y[:,1])
    diff = []
    for i in range(len( eigv )-1):
        diff.append( ( abs(eigv[i] - eigv[i+1]) , eigv[i], eigv[i+1] ) ) 

    diff.sort(key=lambda a: a[0])

    cut_list = []
    for t in diff[-top_k:]:
        C1 = []
        C2 = []
        m = t[1] + (t[2]-t[1])/2
        for n in nodelist:
            index = nodelist.index(n)
            if Y[:,1][index] >= m:
                C1.append(n)
            else:
                C2.append(n)
        cluster_list = [C1, C2]
        C1s = [i for i in start_set if i in C1]
        C1t = [i for i in stop_set if i in C1]
        C2s = [i for i in start_set if i in C2]
        C2t = [i for i in stop_set if i in C2]
        if (C1s !=[] and C1t !=[]) or (C2s !=[] and C2t !=[]):
            continue
        else:
            V = indicator_vector(nodelist, cluster_list)
            c = WCut(D, W, D, V)
            cut_list.append((c,m,C1,C2))

    cut_list.sort(key=lambda a: a[0])
    c, m, C1, C2 = cut_list[0]

    plt.figure(figsize=(15, 8))
    for i in nodelist:
        index = nodelist.index(i)
        if i in C1:
            plt.scatter(x=index, y=Y[:,1][index], color = 'r')
        else:
            plt.scatter(x=index, y=Y[:,1][index], color = 'b')
    plt.title(network_filename + ' second smallest eigenvec with cluster split at ' + str(round(m,4)) + '. WCut='+str(round(c,4)))
    
    if save_file == True:
        plt.savefig(network_filename + 'cluster_split_at_' + str(round(m,4)) + '_cut_'+str(round(c,4))+ '_diagP.png')   # save the figure to file
        plt.close() 
    else:
        plt.show()
    
    return c, m, C1, C2


def create_random_color_pallet(number_of_colors):
    pallet = []
    for j in range(number_of_colors+1):
        rand_colors = ["#"+''.join([random.choice('ABCDEF0123456789') for i in range(6)])]
        pallet.append(rand_colors)
    return pallet

def view_pydot(G, data = None):
    pdot = nx.drawing.nx_pydot.to_pydot(G)

    if data != None:
        for i, edge in enumerate(pdot.get_edges()):
            edge.set_label(round(float(edge.get_attributes()[data]),2))
    plt = Image(pdot.create_png())
    display(plt)

def view_clusters(G, clusters, data = None):
    pallet = create_random_color_pallet(len(clusters))
    color = {}

    for c in clusters:
        for i in G:
            if i in clusters[c]:
                color[i] = pallet[c][0]

    nx.set_node_attributes(G, color, 'color')
    nx.set_node_attributes(G, 'filled', 'style')
    view_pydot(G, data)

def normalized_cut(G, cluster_list, nodelist = None, data = None):

    if nodelist == None:
        nodelist = list(G.nodes())
    W = asym_weight_matrix(G, nodelist, data)
    D = asym_weighted_degree_matrix(G, nodelist, data)

    V = indicator_vector(nodelist, cluster_list)
    cut = WCut(D, W, D, V)
    return cut

def random_cut(N, cluster, nodelist = None, data = None):
    #C1 = random.sample(list(N.nodes()), random.randrange(0, len(N)-1))
    C1 = random.sample(list(N.nodes()), len(cluster))
    C2 = [n for n in list(N.nodes()) if n not in C1]

    cluster_list = [C1, C2]
    c = normalized_cut(N, cluster_list, nodelist, data)

    return c

def random_sample_of_2_way_cut(network_filename, grad_graph_filename,n):
    '''Evaluates the WCut of n random clusters, where the clusters are the same size as the ones calculated
       during the optimized cut.'''

    database_filename = "/home/elizabeth/Desktop/GIT/dsgrn_acdc/networks/" + network_filename + ".db"
    database = Database(database_filename) 
    network_txt_filename = "/home/elizabeth/Desktop/GIT/dsgrn_acdc/networks/" + network_filename + ".txt"

    with open(network_txt_filename,"r") as f:
        network = f.read()

    grad_graph = grad_graph = load_json(grad_graph_filename)

    out_edges = get_number_out_edges_from_string(network)
    FP_Poset = get_FP_Poset(out_edges)[0]
    G = reduce_gradient_graph_to_nodes_of_interest(database, grad_graph, FP_Poset)[0]

    strongcc = strongly_connected_components_by_MGI(G, database)
    cG, scc = condensation(G, strongcc)

    N = nx.DiGraph()
    for node in cG:
        N.add_node(node)
        for edge in cG[node]:
            N.add_edge(node, edge)

    add_source_weight_to_cond(G, N, scc)
    avg_two_way_edge_weight_for_OpenOrd(N, data = 'weight')

    cut, clusters = asym_optimized_normalized_cut(N, 2, nodelist = None, data = 'weight')[1:]
    print('WCut_cG', cut)

    R = [random_cut(N, clusters[0], nodelist = None, data = 'weight') for i in range(0,n)]

    return R