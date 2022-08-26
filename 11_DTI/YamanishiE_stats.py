#Trabajo Fin de Máster - Máster en Métodos Computacionales en ciencias
#Autor: Oier Azurmendi Senar
#Tutores_ Silvestre Vicent Cambra y Mikel Hernaez
#Curso 2021-2022

#Objetivo: Obtener las estadísticas de los grafos creados a partir de la DB de Yamanishi E.

######################################################################################################
import os
import networkx as nx
import pandas as pd 
import numpy as np
from collections import Counter
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = [10, 10]

YE_PATH = '/home/oazurmendis/data/jfuente/DTI/Input4Models/DB/Data/Yamanashi_et_al_GoldStandard/E/interactions/e_admat_dgc_mat_2_line.txt'

# load dtis 
DTI = np.loadtxt(YE_PATH, dtype='str')

# Crear objeto Graph - NetworkX
G = nx.from_edgelist(DTI)
#Cracterísticas del grafo
nx.info(G)

#Diferenciar las drogas y proteínas y obtener las características
drugs, proteins = [], []
for drug, protein in DTI:
    drugs.append(drug), proteins.append(protein)

drugs = list(set(drugs))
proteins = list(set(proteins))
print(f'Info: {nx.info(G)}')
print(f'From those, there are {len(drugs)} drug nodes and {len(proteins)} protein nodes')
print(f'The density of the graph is {nx.density(G):.6f}') # 
print(f'Is directed?: {nx.is_directed(G)}')
# we could select subgraphs with subraph(G, nbunch)

# Obtener Matriz de adyacencia
adj = nx.adjacency_matrix(G)
adj_mat = adj.todense()
adj_mat

#Calculo de la esparsidad
array=adj.toarray()
zeros=np.size(array) - np.count_nonzero(array)
dim=np.size(array)
sparsity=(zeros/dim)*100
print(f'Sparsity percentage: {sparsity}')

#Grafica del grafo de la base de datos
colors = []
shapes = []
for node in G.nodes(data=True):
    if node[0] in drugs:
        nx.set_node_attributes(G, {node[0]: "drug"}, name="ntype")
        colors.append('blue')
        shapes.append('o')
    elif node[0] in proteins:
        nx.set_node_attributes(G, {node[0]: "protein"}, name="ntype")
        colors.append('red')
        shapes.append('s')
    else:
        nx.set_node_attributes(G, {node[0]: "unknown"}, name="ntype")
        colors.append('black')
        shapes.append('X')
        fig, ax = plt.subplots(figsize=(10, 10))

fig, ax = plt.subplots(figsize=(10, 10))
nx.draw(G, 
        with_labels=False, 
        node_color=colors, 
        node_size=15, 
        node_shape='o',
        #pos=layout
       )

font = {"fontname": "Helvetica", "color": "k", "fontweight": "bold", "fontsize": 14}
ax.set_title("Yamanishi E DTI Network", font)
font["color"] = "b"
ax.text(
    0.80,
    0.10,
    "Drugs",
    horizontalalignment="left",
    transform=ax.transAxes,
    fontdict=font,
)
font["color"] = "r"
ax.text(
    0.80,
    0.06,
    "Proteins",
    horizontalalignment="left",
    transform=ax.transAxes,
    fontdict=font,
)
#plt.savefig(f'Yamanishi_E.png')

#Grados del grafo
print(G.degree)

grados = [degree[1] for degree in G.degree]
grados.sort()
print(grados)
Counter(grados)

# Cálculos de las distancias
print(f'Is the graph connected?: {nx.is_connected(G)}')
lcn = list(nx.connected_components(G))
print(f'Number of connected components: {len(lcn)}')

#G.subgraph(lcn[0])
n = 0
#for n in range(len(lcn)):
plt.clf()
subg = G.subgraph(lcn[n])
colors = []
shapes = []
for node in subg.nodes(data=True):
    if node[0] in drugs:
        nx.set_node_attributes(subg, {node[0]: "drug"}, name="ntype")
        colors.append('blue')
        shapes.append('o')
    elif node[0] in proteins:
        nx.set_node_attributes(subg, {node[0]: "protein"}, name="ntype")
        colors.append('red')
        shapes.append('s')
    else:
        nx.set_node_attributes(subg, {node[0]: "unknown"}, name="ntype")
        colors.append('black')
        shapes.append('X')
nx.draw(subg,
        node_size=15,
        node_color=colors)
#plt.savefig(f'YamE/sub_{n}.png')

# Nodos isolados
list(nx.isolates(G))

#Add degree plot
fig, grad = plt.subplots(figsize=(10, 10))
grad.set_title("Degree histogram")
grad.set_xlabel("Degree")
grad.set_ylabel("# of Nodes")
fig.tight_layout()
grad.plot(*np.unique(a, return_counts=True),'o-')
fig.show()
#plt.savefig('Yamanishi-E-Degree-nueva.png')

# for subgraphs
lcn = list(nx.connected_components(G))
size_subgraphs = []
for i in range(len(lcn)):
    subg = G.subgraph(lcn[i])
    size_subgraphs.append(len(subg.nodes))
size_subgraphs.sort()
print(size_subgraphs)

#Nodos por componente
fig, comp = plt.subplots(figsize=(10, 10))
comp.set_title("Connectivity Plot")
comp.set_xlabel("# of Nodes")
comp.set_ylabel("Connected components")
fig.tight_layout()
comp.plot(*np.unique(size_subgraphs, return_counts=True),'o-')
fig.show()
#plt.savefig('Yamanishi-E-Connectivity-nueva.png')
