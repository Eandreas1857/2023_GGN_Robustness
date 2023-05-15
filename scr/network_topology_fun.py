import re
import networkx as nx

def edge_bool_list(string):
    neq_edges = ['~Hb', '~Gt', '~Kr', '~Kni']

    edges = []
    Hb = re.search('Hb :', string).span()
    Gt = re.search('Gt :', string).span()
    Kr = re.search('Kr :', string).span()
    Kni = re.search('Kni :', string).span()
    iter = [('Hb', string[Hb[1]:Gt[0]]), ('Gt', string[Gt[1]:Kr[0]]), ('Kr',string[Kr[1]:Kni[0]]), ('Kni',string[Kni[1]:])]
    for i in iter:
        for ne in neq_edges:
            if ne in i[1]:
                edges.append(((ne[1:], i[0]),0))
            else:
                if ne[1:] in i[1]:
                    edges.append(((ne[1:], i[0]),1))
    return edges

def nxRN(eb_list):
    RN = nx.DiGraph()
    for e, w in eb_list:
        RN.add_edge(e[0], e[1], weight = w)
    return RN

def num_PFL_NFL(RN):    
    #PFL:=cycle with even number of repressing edges
    cycles = sorted(nx.simple_cycles(RN))
    for c in cycles:
        i = c[0]
        c.append(i)

    weight = nx.get_edge_attributes(RN, 'weight')
    PFL = 0 #positive feedback loops
    NFL = 0 #negative feedback loops
    for c in cycles:
        w = 0
        for i in range(len(c)-1):
            w += 1-weight[(c[i], c[i+1])]
        if (w % 2) == 0: #'The number is even'
            PFL += 1
        else: 
            NFL += 1
    return PFL, NFL