#Trabajo Fin de Máster - Máster en Métodos Computacionales en ciencias
#Autor: Oier Azurmendi Senar
#Tutores_ Silvestre Vicent Cambra y Mikel Hernaez
#Curso 2021-2022

#Objetivo: Obtener las gráficas de Connectivity Plot y Degree Plot de todos los DB a la vez para compararlos.

###############################################################################################################
import os
import networkx as nx
import pandas as pd 
import numpy as np
from tdc.multi_pred import DTI
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
plt.rcParams['figure.figsize'] = [10, 10]

#Add paths
BINDINGDB_PATH='/home/oazurmendis/data/jfuente/DTI/Input4Models/DB/Data/BindingDB/BindingDB_All_2021m11/BindingDB_All.tsv'
DRUGBANK_PATH = '/home/oazurmendis/data/jfuente/DTI/Input4Models/DB/Data/DrugBank/DrugBank_DTIs.tsv'
BIOSNAP_PATH = '/home/oazurmendis/data/jfuente/DTI/Input4Models/DB/Data/BIOSNAP/ChG-Miner_miner-chem-gene/ChG-Miner_miner-chem-gene.tsv'
E_PATH = '/home/oazurmendis/data/jfuente/DTI/Input4Models/DB/Data/Yamanashi_et_al_GoldStandard/E/interactions/e_admat_dgc_mat_2_line.txt'
NR_PATH = '/home/oazurmendis/data/jfuente/DTI/Input4Models/DB/Data/Yamanashi_et_al_GoldStandard/NR/interactions/nr_admat_dgc_mat_2_line.txt'
GPCR_PATH = '/home/oazurmendis/data/jfuente/DTI/Input4Models/DB/Data/Yamanashi_et_al_GoldStandard/GPCR/interactions/gpcr_admat_dgc_mat_2_line.txt'
IC_PATH = '/home/oazurmendis/data/jfuente/DTI/Input4Models/DB/Data/Yamanashi_et_al_GoldStandard/IC/interactions/ic_admat_dgc_mat_2_line.txt'

##Cargo las bases de datos

#Davis et al

data = DTI(name = 'DAVIS')
data.binarize(threshold = 30, order = 'descending')
df = data.get_data()
# load dtis 
df = df[df['Y']==1]
DTI_davis = df[['Drug_ID', 'Target_ID']].to_numpy(dtype=str)

# create Graph object - NetworkX
G_davis = nx.from_edgelist(DTI_davis)
DTI_davisdf=pd.DataFrame(DTI_davis)
#DTI_davisdf.to_csv('Davisdf_graph.csv',header=True,index=True)
nx.info(G_davis)

#BindingDB
data = DTI(name = 'BindingDB_Kd')
data.harmonize_affinities(mode = 'max_affinity')
df = data.get_data()
DTI_bindingdb = df[['Drug_ID', 'Target_ID']].to_numpy(dtype=str)
G_bindingdb = nx.from_edgelist(DTI_bindingdb)
DTI_bindingdbdf=pd.DataFrame(DTI_bindingdb)
#DTI_bindingdbdf.to_csv('Bindingdbdf_graph.csv',header=True,index=True)
nx.info(G_bindingdb)

# Drugbank
DTI_drugbank = np.loadtxt(DRUGBANK_PATH, dtype='str')
# create Graph object - NetworkX
G_drugbank = nx.from_edgelist(DTI_drugbank)
DTI_drugbankdf=pd.DataFrame(DTI_drugbank)
#DTI_drugbankdf.to_csv('Drugbank_graph.csv',header=True,index=True)
nx.info(G_drugbank)

# BIOSNAP
DTI_biosnap = np.loadtxt(BIOSNAP_PATH, dtype='str')
# create Graph object - NetworkX
G_biosnap = nx.from_edgelist(DTI_biosnap)
DTI_biosnapdf=pd.DataFrame(DTI_biosnap)
#DTI_biosnapdf.to_csv('Biosnap_graph.csv',header=True,index=True)
nx.info(G_biosnap)

# Yamanishi E
DTI_e = np.loadtxt(E_PATH, dtype='str')
# create Graph object - NetworkX
G_e = nx.from_edgelist(DTI_e)
DTI_edf=pd.DataFrame(DTI_e)
#DTI_edf.to_csv('Edf_graph.csv',header=True,index=True)
nx.info(G_e)

# Yamanishi NR
DTI_nr = np.loadtxt(NR_PATH, dtype='str')
# create Graph object - NetworkX
G_nr = nx.from_edgelist(DTI_nr)
DTI_nrdf=pd.DataFrame(DTI_nr)
#DTI_nrdf.to_csv('NRdf_graph.csv',header=True,index=True)
nx.info(G_nr)

# Yamanishi GPCR
DTI_gpcr = np.loadtxt(GPCR_PATH, dtype='str')
# create Graph object - NetworkX
G_gpcr = nx.from_edgelist(DTI_gpcr)
DTI_gpcrdf=pd.DataFrame(DTI_gpcr)
#DTI_gpcrdf.to_csv('GPCRdf_graph.csv',header=True,index=True)
nx.info(G_gpcr)

# Yamanishi IC
DTI_ic = np.loadtxt(IC_PATH, dtype='str')
# create Graph object - NetworkX
G_ic = nx.from_edgelist(DTI_ic)
DTI_icdf=pd.DataFrame(DTI_ic)
#DTI_icdf.to_csv('ICdf_graph.csv',header=True,index=True)
nx.info(G_ic)

##Obtención de los degrees de cada DB

#Davis et al
degree_davis = [degree[1] for degree in G_davis.degree]
degree_davis.sort()
Davis_Dregree=pd.DataFrame(np.unique(degree_davis))
#Davis_Dregree.to_csv('Davis_Dregree_graph.csv',header=True,index=True)

#BindingDB
degree_bindingdb=[degree[1] for degree in G_bindingdb.degree]
degree_bindingdb.sort()
BindingDB_Dregree=pd.DataFrame(np.unique(degree_bindingdb))
#BindingDB_Dregree.to_csv('BindingDB_Dregree_graph.csv',header=True,index=True)

#Biosnap
degree_biosnap=[degree[1] for degree in G_biosnap.degree]
degree_biosnap.sort()
biosnap_Dregree=pd.DataFrame(np.unique(degree_biosnap))
#biosnap_Dregree.to_csv('Biosnap_Dregree_graph.csv',header=True,index=True)

#Drugbank
degree_drugbank=[degree[1] for degree in G_drugbank.degree]
degree_drugbank.sort()
Drugbank_Dregree=pd.DataFrame(np.unique(degree_drugbank))
#Drugbank_Dregree.to_csv('Drugbank_Dregree_graph.csv',header=True,index=True)

#Yamanishi E
degree_e=[degree[1] for degree in G_e.degree]
degree_e.sort()
E_Dregree=pd.DataFrame(np.unique(degree_e))
#E_Dregree.to_csv('E_Dregree_graph.csv',header=True,index=True)

#Yamanishi NR
degree_nr=[degree[1] for degree in G_nr.degree]
degree_nr.sort()
NR_Dregree=pd.DataFrame(np.unique(degree_nr))
#NR_Dregree.to_csv('NR_Dregree_graph.csv',header=True,index=True)

#Yamanishi GPCR
degree_gpcr=[degree[1] for degree in G_gpcr.degree]
degree_gpcr.sort()
GPCR_Dregree=pd.DataFrame(np.unique(degree_gpcr))
#GPCR_Dregree.to_csv('GPCR_Dregree_graph.csv',header=True,index=True)

#Yamanishi IC
degree_ic=[degree[1] for degree in G_ic.degree]
degree_ic.sort()
IC_Dregree=pd.DataFrame(np.unique(degree_ic))
#IC_Dregree.to_csv('IC_Dregree_graph.csv',header=True,index=True)

## Degree Plots
#Se intenta ver qué escala representa mejor de una forma visual todos los DB.

#Escala logaritmica base 2
fig, ax = plt.subplots(figsize=(10, 10))
ax.set_title("Degree histogram in log")
ax.set_xlabel("Degree")
ax.set_ylabel("# of Nodes")
fig.tight_layout()
ax.plot(*np.unique(np.log(degree_davis), return_counts=True),'o-',label='Davis')
ax.plot(*np.unique(np.log(degree_bindingdb), return_counts=True),'o-',label='BindingDB')
ax.plot(*np.unique(np.log(degree_biosnap), return_counts=True),'o-',label='Biosnap')
ax.plot(*np.unique(np.log(degree_drugbank), return_counts=True),'o-',label='Drugbank')
ax.plot(*np.unique(np.log(degree_e), return_counts=True),'o-',label='Yamanishi E')
ax.plot(*np.unique(np.log(degree_nr), return_counts=True),'o-',label='Yamanishi NR')
ax.plot(*np.unique(np.log(degree_gpcr), return_counts=True),'o-',label='Yamanishi GPCR')
ax.plot(*np.unique(np.log(degree_ic), return_counts=True),'o-',label='Yamanishi IC')
ax.legend()
plt.show()
#ax.set_xlim(0,25)
#plt.savefig('Degree-log.png')

#Escala semi-logarítmica
fig, ax = plt.subplots(figsize=(10, 10))
ax.set_title("Degree histogram in log (Y-axis)")
ax.set_xlabel("Degree")
ax.set_ylabel("Log of # of Nodes")
fig.tight_layout()
ax.semilogy(*np.unique(degree_davis, return_counts=True),'o-',label='Davis')
ax.semilogy(*np.unique(degree_bindingdb, return_counts=True),'o-',label='BindingDB')
ax.semilogy(*np.unique(degree_biosnap, return_counts=True),'o-',label='Biosnap')
ax.semilogy(*np.unique(degree_drugbank, return_counts=True),'o-',label='Drugbank')
ax.semilogy(*np.unique(degree_e, return_counts=True),'o-',label='Yamanishi E')
ax.semilogy(*np.unique(degree_nr, return_counts=True),'o-',label='Yamanishi NR')
ax.semilogy(*np.unique(degree_gpcr, return_counts=True),'o-',label='Yamanishi GPCR')
ax.semilogy(*np.unique(degree_ic, return_counts=True),'o-',label='Yamanishi IC')
ax.legend()
plt.show()
#plt.savefig('Degree.png')
#ax.set_xlim(0,25)
#plt.savefig('Degree-semilog.png')

#Escala logaritmica base 10
fig, ax = plt.subplots(figsize=(10, 10))
ax.set_title("Degree histogram in log")
ax.set_xlabel("Degree")
ax.set_ylabel("# of Nodes")
fig.tight_layout()
ax.loglog(*np.unique(degree_davis, return_counts=True),'o-',label='Davis')
ax.loglog(*np.unique(degree_bindingdb, return_counts=True),'o-',label='BindingDB')
ax.loglog(*np.unique(degree_biosnap, return_counts=True),'o-',label='Biosnap')
ax.loglog(*np.unique(degree_drugbank, return_counts=True),'o-',label='Drugbank')
ax.loglog(*np.unique(degree_e, return_counts=True),'o-',label='Yamanishi E')
ax.loglog(*np.unique(degree_nr, return_counts=True),'o-',label='Yamanishi NR')
ax.loglog(*np.unique(degree_gpcr, return_counts=True),'o-',label='Yamanishi GPCR')
ax.loglog(*np.unique(degree_ic, return_counts=True),'o-',label='Yamanishi IC')
ax.legend()
plt.show()
#plt.savefig('Degree.png')
#ax.set_xlim(0,25)
#plt.savefig('Degree-log_10scale.png')

# Aplicando el zoom en la región interesada
fig, ax = plt.subplots(figsize=(10, 10))
axins = inset_axes(ax, 3.5,4 , loc=9)
ax.set_title("Degree histogram in log")
ax.set_xlabel("Degree")
ax.set_ylabel("# of Nodes")
fig.tight_layout()
ax.loglog(*np.unique(degree_davis, return_counts=True),'o-',label='Davis')
ax.loglog(*np.unique(degree_bindingdb, return_counts=True),'o-',label='BindingDB')
ax.loglog(*np.unique(degree_biosnap, return_counts=True),'o-',label='Biosnap')
ax.loglog(*np.unique(degree_drugbank, return_counts=True),'o-',label='Drugbank')
ax.loglog(*np.unique(degree_e, return_counts=True),'o-',label='Yamanishi E')
ax.loglog(*np.unique(degree_nr, return_counts=True),'o-',label='Yamanishi NR')
ax.loglog(*np.unique(degree_gpcr, return_counts=True),'o-',label='Yamanishi GPCR')
ax.loglog(*np.unique(degree_ic, return_counts=True),'o-',label='Yamanishi IC')
axins.loglog(*np.unique(degree_davis, return_counts=True),'o-',label='Davis')
axins.loglog(*np.unique(degree_bindingdb, return_counts=True),'o-',label='BindingDB')
axins.loglog(*np.unique(degree_biosnap, return_counts=True),'o-',label='Biosnap')
axins.loglog(*np.unique(degree_drugbank, return_counts=True),'o-',label='Drugbank')
axins.loglog(*np.unique(degree_e, return_counts=True),'o-',label='Yamanishi E')
axins.loglog(*np.unique(degree_nr, return_counts=True),'o-',label='Yamanishi NR')
axins.loglog(*np.unique(degree_gpcr, return_counts=True),'o-',label='Yamanishi GPCR')
axins.loglog(*np.unique(degree_ic, return_counts=True),'o-',label='Yamanishi IC')
axins.set_xlim(10, 100) # apply the x-limits
axins.set_ylim(0, 150)
mark_inset(ax, axins, loc1=3, loc2=4, fc="grey", ec="0.5")
ax.legend()
plt.show()
#plt.savefig('Degree-log_zoom.png')


##Connectivity 

#Obtienen valores de los componentes conectados en cada DB

# Davis et al
lcn = list(nx.connected_components(G_davis))
size_subgraphs_davis = []
for i in range(len(lcn)):
    subg = G_davis.subgraph(lcn[i])
    size_subgraphs_davis.append(len(subg.nodes))
size_subgraphs_davis.sort()
Davis_Subgraph=pd.DataFrame(size_subgraphs_davis)
#Davis_Subgraph.to_csv('Connectivity_davis.csv',header=True,index=True)

#BindingDB
lcn = list(nx.connected_components(G_bindingdb))
size_subgraphs_bindingdb = []
for i in range(len(lcn)):
    subg = G_bindingdb.subgraph(lcn[i])
    size_subgraphs_bindingdb.append(len(subg.nodes))
size_subgraphs_bindingdb.sort()
bindingdb_Subgraph=pd.DataFrame(size_subgraphs_bindingdb)
#bindingdb_Subgraph.to_csv('Connectivity_bindingdb.csv',header=True,index=True)

#Biosnap
lcn = list(nx.connected_components(G_biosnap))
size_subgraphs_biosnap = []
for i in range(len(lcn)):
    subg = G_biosnap.subgraph(lcn[i])
    size_subgraphs_biosnap.append(len(subg.nodes))
size_subgraphs_biosnap.sort()
biosnap_Subgraph=pd.DataFrame(size_subgraphs_biosnap)
#biosnap_Subgraph.to_csv('Connectivity_biosnap.csv',header=True,index=True)

#Drugbank
lcn = list(nx.connected_components(G_drugbank))
size_subgraphs_drugbank= []
for i in range(len(lcn)):
    subg = G_drugbank.subgraph(lcn[i])
    size_subgraphs_drugbank.append(len(subg.nodes))
size_subgraphs_drugbank.sort()
drugbank_Subgraph=pd.DataFrame(size_subgraphs_drugbank)
#drugbank_Subgraph.to_csv('Connectivity_drugbank.csv',header=True,index=True)

#Yamanishi E
lcn = list(nx.connected_components(G_e))
size_subgraphs_e = []
for i in range(len(lcn)):
    subg = G_e.subgraph(lcn[i])
    size_subgraphs_e.append(len(subg.nodes))
size_subgraphs_e.sort()
e_Subgraph=pd.DataFrame(size_subgraphs_e)
#e_Subgraph.to_csv('Connectivity_e.csv',header=True,index=True)

#Yamanishi GR
lcn = list(nx.connected_components(G_nr))
size_subgraphs_nr = []
for i in range(len(lcn)):
    subg = G_nr.subgraph(lcn[i])
    size_subgraphs_nr.append(len(subg.nodes))
size_subgraphs_nr.sort()
nr_Subgraph=pd.DataFrame(size_subgraphs_nr)
#nr_Subgraph.to_csv('Connectivity_nr.csv',header=True,index=True)

#Yamanishi GPCR
lcn = list(nx.connected_components(G_gpcr))
size_subgraphs_gpcr = []
for i in range(len(lcn)):
    subg = G_gpcr.subgraph(lcn[i])
    size_subgraphs_gpcr.append(len(subg.nodes))
size_subgraphs_gpcr.sort()
gpcr_Subgraph=pd.DataFrame(size_subgraphs_gpcr)
#gpcr_Subgraph.to_csv('Connectivity_gpcr.csv',header=True,index=True)

#Yamanishi IC
lcn = list(nx.connected_components(G_ic))
size_subgraphs_ic = []
for i in range(len(lcn)):
    subg = G_ic.subgraph(lcn[i])
    size_subgraphs_ic.append(len(subg.nodes))
size_subgraphs_ic.sort()
ic_Subgraph=pd.DataFrame(size_subgraphs_ic)
#ic_Subgraph.to_csv('Connectivity_ic.csv',header=True,index=True)

##Connectivity plots
#Ver la cantidad de nodos por componente conectado
#Se intenta ver que escala se le ajusta mejor

fig, ax3 = plt.subplots(figsize=(10, 10))
ax3.set_title("Connectivity Plot in Log")
ax3.set_xlabel("# of Nodes")
ax3.set_ylabel("Connected components")
fig.tight_layout()
ax3.plot(*np.unique(np.log(size_subgraphs_davis), return_counts=True),'o-',label='Davis')
ax3.plot(*np.unique(np.log(size_subgraphs_bindingdb), return_counts=True),'o-',label='BindingDB')
ax3.plot(*np.unique(np.log(size_subgraphs_biosnap), return_counts=True),'o-',label='Biosnap')
ax3.plot(*np.unique(np.log(size_subgraphs_drugbank), return_counts=True),'o-',label='Drugbank')
ax3.plot(*np.unique(np.log(size_subgraphs_e), return_counts=True),'o-',label='Yamanishi E')
ax3.plot(*np.unique(np.log(size_subgraphs_nr), return_counts=True),'o-',label='Yamanishi NR')
ax3.plot(*np.unique(np.log(size_subgraphs_gpcr), return_counts=True),'o-',label='Yamanishi GPCR')
ax3.plot(*np.unique(np.log(size_subgraphs_ic), return_counts=True),'o-',label='Yamanishi IC')
ax3.legend()
fig.show()
#plt.savefig('Connectivity.png')
#ax3.set_xlim(0,25)
#plt.savefig('Connectivity-log.png')

#Escala semilogarítmica
fig, ax3 = plt.subplots(figsize=(10, 10))
ax3.set_title("Connectivity Plot in Log (Y-axis)")
ax3.set_xlabel("# of Nodes")
ax3.set_ylabel("Log of Connected components")
fig.tight_layout()
ax3.semilogy(*np.unique(size_subgraphs_davis, return_counts=True),'o-',label='Davis')
ax3.semilogy(*np.unique(size_subgraphs_bindingdb, return_counts=True),'o-',label='BindingDB')
ax3.semilogy(*np.unique(size_subgraphs_biosnap, return_counts=True),'o-',label='Biosnap')
ax3.semilogy(*np.unique(size_subgraphs_drugbank, return_counts=True),'o-',label='Drugbank')
ax3.semilogy(*np.unique(size_subgraphs_e, return_counts=True),'o-',label='Yamanishi E')
ax3.semilogy(*np.unique(size_subgraphs_nr, return_counts=True),'o-',label='Yamanishi NR')
ax3.semilogy(*np.unique(size_subgraphs_gpcr, return_counts=True),'o-',label='Yamanishi GPCR')
ax3.semilogy(*np.unique(size_subgraphs_ic, return_counts=True),'o-',label='Yamanishi IC')
ax3.legend()
fig.show()
#plt.savefig('Connectivity.png')
#ax3.set_xlim(0,25)
#plt.savefig('Connectivity-semilog.png')

#Escala logaritmica en base 10
fig, ax3 = plt.subplots(figsize=(10, 10))
ax3.set_title("Connectivity Plot in Log")
ax3.set_xlabel("# of Nodes")
ax3.set_ylabel("Connected components")
fig.tight_layout()
ax3.loglog(*np.unique(size_subgraphs_davis, return_counts=True),'o-',label='Davis')
ax3.loglog(*np.unique(size_subgraphs_bindingdb, return_counts=True),'o-',label='BindingDB')
ax3.loglog(*np.unique(size_subgraphs_biosnap, return_counts=True),'o-',label='Biosnap')
ax3.loglog(*np.unique(size_subgraphs_drugbank, return_counts=True),'o-',label='Drugbank')
ax3.loglog(*np.unique(size_subgraphs_e, return_counts=True),'o-',label='Yamanishi E')
ax3.loglog(*np.unique(size_subgraphs_nr, return_counts=True),'o-',label='Yamanishi NR')
ax3.loglog(*np.unique(size_subgraphs_gpcr, return_counts=True),'o-',label='Yamanishi GPCR')
ax3.loglog(*np.unique(size_subgraphs_ic, return_counts=True),'o-',label='Yamanishi IC')
ax3.legend()
fig.show()
#plt.savefig('Connectivity.png')
#plt.savefig('Connectivity-log_10_scale.png')

#Escala logaritmica en base 10 con el zoom activado
fig, ax3 = plt.subplots(figsize=(10, 10))
axins1 = inset_axes(ax3, 4,4 , loc=9)
ax3.set_title("Connectivity Plot in Log")
ax3.set_xlabel("# of Nodes")
ax3.set_ylabel("Connected components")
fig.tight_layout()
ax3.loglog(*np.unique(size_subgraphs_davis, return_counts=True),'o-',label='Davis')
ax3.loglog(*np.unique(size_subgraphs_bindingdb, return_counts=True),'o-',label='BindingDB')
ax3.loglog(*np.unique(size_subgraphs_biosnap, return_counts=True),'o-',label='Biosnap')
ax3.loglog(*np.unique(size_subgraphs_drugbank, return_counts=True),'o-',label='Drugbank')
ax3.loglog(*np.unique(size_subgraphs_e, return_counts=True),'o-',label='Yamanishi E')
ax3.loglog(*np.unique(size_subgraphs_nr, return_counts=True),'o-',label='Yamanishi NR')
ax3.loglog(*np.unique(size_subgraphs_gpcr, return_counts=True),'o-',label='Yamanishi GPCR')
ax3.loglog(*np.unique(size_subgraphs_ic, return_counts=True),'o-',label='Yamanishi IC')
axins1.loglog(*np.unique(size_subgraphs_davis, return_counts=True),'o-',label='Davis')
axins1.loglog(*np.unique(size_subgraphs_bindingdb, return_counts=True),'o-',label='BindingDB')
axins1.loglog(*np.unique(size_subgraphs_biosnap, return_counts=True),'o-',label='Biosnap')
axins1.loglog(*np.unique(size_subgraphs_drugbank, return_counts=True),'o-',label='Drugbank')
axins1.loglog(*np.unique(size_subgraphs_e, return_counts=True),'o-',label='Yamanishi E')
axins1.loglog(*np.unique(size_subgraphs_nr, return_counts=True),'o-',label='Yamanishi NR')
axins1.loglog(*np.unique(size_subgraphs_gpcr, return_counts=True),'o-',label='Yamanishi GPCR')
axins1.loglog(*np.unique(size_subgraphs_ic, return_counts=True),'o-',label='Yamanishi IC')
axins1.set_xlim(3, 30)
axins1.set_ylim(0, 75) 
mark_inset(ax3, axins1, loc1=3, loc2=4, fc="grey", ec="0.5")
ax3.legend()
fig.show()
#plt.savefig('Connectivity-zoom.png')

