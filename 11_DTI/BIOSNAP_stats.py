#Trabajo Fin de Máster - Máster en Métodos Computacionales en ciencias
#Autor: Oier Azurmendi Senar
#Tutores_ Silvestre Vicent Cambra y Mikel Hernaez
#Curso 2021-2022

#Objetivo: Obtener las estadísticas de los grafos creados a partir de la DB de Biosnap.

######################################################################################################
import os
import networkx as nx
import pandas as pd 
import numpy as np
from collections import Counter
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = [10, 10]

BIOSNAP_PATH = '/home/oazurmendis/data/jfuente/DTI/Input4Models/DB/Data/BIOSNAP/ChG-Miner_miner-chem-gene/ChG-Miner_miner-chem-gene.tsv'

# load dtis 
DTI = np.loadtxt(BIOSNAP_PATH, dtype='str')
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

#Cálculo de la esparsidad
array=adj.toarray()
zeros=np.size(array) - np.count_nonzero(array)
dim=np.size(array)
sparsity=(zeros/dim)*100
print(f'Sparsity percentage: {sparsity}')

#Grafica del grafo de la base de datos
fig, ax = plt.subplots(figsize=(10, 10))

nx.draw(G, 
        with_labels=False, 
        node_color=colors, 
        node_size=15, 
        node_shape='o',
        #pos=layout
       )

font = {"fontname": "Helvetica", "color": "k", "fontweight": "bold", "fontsize": 14}
ax.set_title("BIOSNAP DTI Network", font)
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
#plt.savefig(f'BIOSNAP.png')

#Grados del grafo
print(G.degree)

grados = [degree[1] for degree in G.degree]
grados.sort()
print(grados)
Counter(grados)

# Cálculos de la distancia
print(f'Is the graph connected?: {nx.is_connected(G)}')
lcn = list(nx.connected_components(G))
print(f'Number of connected components: {len(lcn)}')

#Componentes conectados
fig, comp = plt.subplots(figsize=(10, 10))
Gcc = G.subgraph(sorted(nx.connected_components(G), key=len, reverse=True)[0])
pos = nx.spring_layout(Gcc, seed=10396953)
nx.draw_networkx_nodes(Gcc, pos, ax=ax0, node_size=20)
nx.draw_networkx_edges(Gcc, pos, ax=ax0, alpha=0.4)
comp.set_title("Connected components of G")
fig.tight_layout()
plt.show()

fig, grad = plt.subplots(figsize=(10, 10))
grad.plot(grados, "b-", marker="o")
grad.set_title("Degree Rank Plot")
grad.set_ylabel("Degree")
grad.set_xlabel("Rank")
fig.tight_layout()
plt.show()

fig, ax2 = plt.subplots(figsize=(10, 10))
ax2.set_title("Degree histogram")
ax2.set_xlabel("Degree")
ax2.set_ylabel("# of Nodes")
fig.tight_layout()
ax2.plot(*np.unique(grados, return_counts=True),'o-')
plt.show()
plt.savefig('BIOSNAP-Degree-nueva.png')

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
#plt.savefig(f'Biosnap/sub_{n}.png')

# Nodos izolados
print(list(nx.isolates(G)))

# for subgraphs
lcn = list(nx.connected_components(G))
size_subgraphs = []
for i in range(len(lcn)):
    subg = G.subgraph(lcn[i])
    size_subgraphs.append(len(subg.nodes))

size_subgraphs.sort()
print(size_subgraphs)

#Nodos por componente
fig, ax3 = plt.subplots(figsize=(10, 10))
ax3.set_title("Connectivity Plot")
ax3.set_xlabel("# of Nodes")
ax3.set_ylabel("Connected components")
fig.tight_layout()
ax3.plot(*np.unique(size_subgraphs, return_counts=True),'o-')
fig.show()
plt.savefig('BIOSNAP-Connectivity-nueva.png')
