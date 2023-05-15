import networkx as nx

def score_region_transitions(database, G, FP_Regions, scc = None):
    '''
    Calculating a 'score' for each region transition by counting the number of edges from region r to r+1 
    and norming it by the number of edges out of region r total. If each edge in the gradient graph has the same 
    chance of being take, then we can think of this score as the probability of transitioning from region r to r+1.
    returns: average 'score' from each region transition. 1 is best possible score, 0 is the worst. Also return a
    probabily of taking path r = 1, 2, 3, 4, 5, 6, 7, 8. Additionally, it scores all paths from region 1 to 8 and 
    returns the percentile that the path 1, 2, 3, 4, 5, 6, 7, 8 falls in when the paths are ranked by score. 1.0 
    means it has the highest score, 0 means it has the lowest.
    '''
    R = {}
    for i in FP_Regions:
        for j in FP_Regions[i]:
            R[j] = i

    c = database.conn.cursor()

    Region_count = {}
    for edge in G.edges():
        if scc == None:
            s = edge[0][-1]
            t = edge[1][-1]
        else:
            s = scc[edge[0]][0][-1]
            t = scc[edge[1]][0][-1]

        sMGI = c.execute('select MorseGraphIndex from Signatures where ParameterIndex is ' + str(s))
        MGI = sMGI.fetchone()[0]
        sFP = [row[0] for row in c.execute('select Label from MorseGraphAnnotations where MorseGraphIndex is ' + str(MGI))]

        tMGI = c.execute('select MorseGraphIndex from Signatures where ParameterIndex is ' + str(t))
        MGI = tMGI.fetchone()[0]
        tFP = [row[0] for row in c.execute('select Label from MorseGraphAnnotations where MorseGraphIndex is ' + str(MGI))]

        if R[sFP[0]]!=R[tFP[0]]:
            if (R[sFP[0]], R[tFP[0]]) not in Region_count:
                Region_count[(R[sFP[0]], R[tFP[0]])] = 1
            else:
                Region_count[(R[sFP[0]], R[tFP[0]])] += 1


    out = {}
    for i in range(1,9):
        out[i] = 0
        for r in Region_count:
            if r[0] == i:
                out[i] += Region_count[r]
    avg = 0
    prob = 1
    for i in range(1,8):
        try:
            print((i,i+1), Region_count[(i,i+1)]/out[i])
            prob = prob*(Region_count[(i,i+1)]/out[i])
            avg += Region_count[(i,i+1)]/out[i]
        except:
            continue
    RG = nx.DiGraph()

    for i in Region_count:
        RG.add_edge(i[0], i[1], weight = Region_count[i]/out[i[0]]) 

    simple = [p for p in nx.all_simple_paths(RG, 1, 8)]

    tup = []
    for p in simple:
        score = 0
        for n in range(len(p)-1):
            score += Region_count[(p[n], p[n+1])]/out[p[n]]
        tup.append((p,score/7))
    tup.sort(key = lambda x: x[1])
    for t in tup:
        if t[0] == [1, 2, 3, 4, 5, 6, 7, 8]:
            percentile = tup.index(t)/(len(tup)-1)
        
    return avg/7, prob, percentile