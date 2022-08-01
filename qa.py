import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import string
import urllib
import tsplib95


url = "https://people.sc.fsu.edu/~jburkardt/datasets/tsp/gr17.tsp"
file = urllib.request.urlopen(url)
tsp_graph = file.read().decode('utf-8')
tsp_graph=tsplib95.parse(tsp_graph)
graph = tsp_graph.get_graph() 
mapping = dict(zip(graph, string.ascii_uppercase))
graph = nx.relabel_nodes(graph, mapping)
edges = list(graph.edges)
nodes = list(graph.nodes)
cmap = ['blue' for node in nodes]
pos = nx.spring_layout(graph, seed=7)
nx.draw(graph, pos, with_labels=True)
labels = nx.get_edge_attributes(graph,'weight')

nx.draw_networkx_edge_labels(graph,pos,edge_labels=labels)



from collections import defaultdict
from itertools import permutations
h = defaultdict(float)
J = defaultdict(float)
C = graph.size(weight="weight")*graph.number_of_nodes()/graph.number_of_edges()
#C = 1000000000
N = len(nodes)

#nodes constraint
for node in nodes:
    for j in range(N):
        J[(node, j),(node, j)] -= C
        for j_1 in range(j+1,N):
            J[(node, j), (node, j_1)] += 2*C
              
#order constraint
for j in range(N):
    for i, node in enumerate(nodes):
        J[(node, j),(node, j)] -= C
        for node_1 in nodes[i+1:]:
            J[(node, j), (node_1, j)] += 2*C

#connection constraint
comb = permutations(nodes, 2)
unvalid_conn = set(edges) - set(comb)
for i in unvalid_conn:
    for j in range(N):
        J[(i[0], j), (i[1], j+1%N)] += C*0.5
        
#cost funtion
for edge in edges:
    weight = graph.get_edge_data(*edge)['weight']
    for j in range(N):
        J[(edge[0], j), (edge[1], (j+1)%N)] += weight
        J[(edge[1], j), (edge[0], (j+1)%N)] += weight



from dimod import BinaryQuadraticModel
from dwave.system import EmbeddingComposite, DWaveSampler
from hybrid.reference.kerberos import KerberosSampler
sampler = EmbeddingComposite(DWaveSampler())
bqm = BinaryQuadraticModel(J, offset = 2*C*N, vartype = 'BINARY')
sampleset = KerberosSampler().sample(bqm, num_reads = 10)


sample = sampleset.first.sample

route = [None]*(N+1)
for el in sample.items():
        #print(el)
        if el[1]:
            route[el[0][1]] = el[0][0]

path = [None]*N
for t in range(N):
    path[t] = (route[t], route[(t+1)%N], graph.edges[route[t], route[(t+1)%N]]['weight'])
path

tsp_graph = nx.Graph()
tsp_graph.add_weighted_edges_from(path)
pos = nx.spring_layout(graph, seed=7)
nx.draw(tsp_graph, pos, with_labels=True)
labels = nx.get_edge_attributes(tsp_graph,'weight')
nx.draw_networkx_edge_labels(tsp_graph,pos,edge_labels=labels)

tsp_graph.size('weight')