#Trabajo Fin de Máster - Máster en Métodos Computacionales en ciencias
#Autor: Oier Azurmendi Senar
#Tutores_ Silvestre Vicent Cambra y Mikel Hernaez
#Curso 2021-2022

#Objetivo: Obtener las estadísticas de los grafos creados a partir de la DB de Drugbank.

######################################################################################################
import os
import networkx as nx
import pandas as pd 
import numpy as np
from collections import Counter
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = [10, 10]

Drugbank_PATH = '/home/oazurmendis/data/jfuente/DTI/Input4Models/DB/Data/DrugBank/DrugBank_DTIs.tsv'
# load dtis 
DTI = np.loadtxt(Drugbank__PATH, dtype='str')
# Crear objeto Graph - NetworkX
G = nx.from_edgelist(DTI)
#Caracteristicas del grafo
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

#Cálculo de la esparsidad
array=adj.toarray()
zeros=np.size(array) - np.count_nonzero(array)
dim=np.size(array)
sparsity=(zeros/dim)*100
print(f'Sparsity percentage: {sparsity}')


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

#Gráfica del grado de la base de datos        
fig, ax = plt.subplots(figsize=(10, 10))
nx.draw(G, 
        with_labels=False, 
        node_color=colors, 
        node_size=15, 
        node_shape='o',
        #pos=layout
       )

font = {"fontname": "Helvetica", "color": "k", "fontweight": "bold", "fontsize": 14}
ax.set_title("Drugbank DTI Network", font)
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
#plt.savefig(f'Drugbank.png')

#Grados del grafo
print(G.degree)

grados = [degree[1] for degree in G.degree]
grados.sort()
print(grados)
Counter(grados)

#Cálculos de la distancia
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
#plt.savefig(f'Drugbank/sub_{n}.png')

# isolated nodes
list(nx.isolates(G))

#Gráfica de los grados

fig, grad = plt.subplots(figsize=(10, 10))
grad.set_title("Degree histogram")
grad.set_xlabel("Degree")
grad.set_ylabel("# of Nodes")
fig.tight_layout()
grad.plot(*np.unique(a, return_counts=True),'o-')
fig.show()
plt.savefig('Drugbank-Degree-nueva.png')


# for subgraphs
lcn = list(nx.connected_components(G))
size_subgraphs = []
for i in range(len(lcn)):
    subg = G.subgraph(lcn[i])
    size_subgraphs.append(len(subg.nodes))
size_subgraphs.sort()
print(size_subgraphs)

#Nodos por componente
fig, componentes = plt.subplots(figsize=(10, 10))
componentes.set_title("Connectivity Plot")
componentes.set_xlabel("# of Nodes")
componentes.set_ylabel("Connected components")
fig.tight_layout()
componentes.plot(*np.unique(size_subgraphs, return_counts=True),'o-')
fig.show()
#plt.savefig('Drugbank-Connectivity-nueva.png')

