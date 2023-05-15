import networkx as nx

def create_grad_graph_w_subgraphs_graphml(database, scc, cG, path_graph, prod_graph, FP_Region, start_set, stop_set, Filename):
    ''' graphml filetype '''
    c = database.conn.cursor()
    
    Kni_att = {} # FG layer for Kni node belongs too
    Hb_att = {} # FG layer for Hb node belongs too
    MGI_att = {} # what is the Morse Graph index?
    Region_att = {} # what region if FP associated with 
    
    graph = {} # condenced chemical gradient graph only or apart of path graph? 
    s_t = {} # set of starting and stopping nodes


    for node in cG:
        p = scc[node][0][-1]
        MGI_result = c.execute('select MorseGraphIndex from Signatures where ParameterIndex is ' + str(p))
        MGI = MGI_result.fetchone()[0]

        MGI_att[node] = MGI
        FP_result = [row[0] for row in c.execute('select Label from MorseGraphAnnotations where MorseGraphIndex is ' + str(MGI))]
        if len(FP_result) == 1:
            for r in FP_Region:
                if FP_result[0] in FP_Region[r]:
                    Region_att[node] = r
        else:
            Region_att[node] = 'not mono-stable'
        
        Hb_att[node] = scc[node][0][0]
        Kni_att[node] = scc[node][0][1]

        if node in path_graph:
            graph[node] = 'path'
        elif node in prod_graph:
            graph[node] = 'product'
        elif node in cG:
            graph[node] = 'cG'

    for node in start_set:
        s_t[node] = 'starting'
    for node in stop_set:
        s_t[node] = 'stopping'

    nx.set_node_attributes(cG, Hb_att, 'Hb_FG_layer')
    nx.set_node_attributes(cG, Kni_att, 'Kni_FG_layer')
    nx.set_node_attributes(cG, MGI_att, 'MGI')
    nx.set_node_attributes(cG, Region_att, 'Region')
    nx.set_node_attributes(cG, graph, 'group')
    
    nx.set_node_attributes(cG, s_t, 'start_stop')

    nx.write_graphml(cG, Filename)
