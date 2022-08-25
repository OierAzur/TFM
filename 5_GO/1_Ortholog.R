#Trabajo Fin de Máster - Máster en Métodos Computacionales en ciencias
#Autor: Oier Azurmendi Senar
#Tutores_ Silvestre Vicent Cambra y Mikel Hernaez
#Curso 2021-2022

#Objetivo: Obtener los ortógos de humano para poder realizar los análisis de Gene Ontology

######################################################################################################
setwd("Documents/TFM/GO/")
library(babelgene)
load("/Users/oierazur/TFM/DGE/CCARag_vs_normal_filtro1.Rdata")
load("/Users/oierazur/TFM/DGE/Dieta_vs_normal_filtro1.Rdata")
load("/Users/oierazur/TFM/DGE/Hidrodinamica_vs_normalfiltro1.Rdata")
load("/Users/oierazur/TFM/DGE/Transformadas_vs_normalfiltro1.Rdata")
load("/Users/oierazur/TFM/DGE/Dieta_vs_Hidrodinamica.Rdata")

#CCA (Raf-F1)
ortolog_CCARag <- orthologs(genes=CCARag_vs_normal$ID,species="Mus musculus",human=F)
CCARag_vs_normal_ortolog <- rename(CCARag_vs_normal,symbol="ID")
CCARag_vs_normal_ortolog <- merge(ortolog_CCARag,CCARag_vs_normal_ortolog,
                                  by="symbol") 
save(CCARag_vs_normal_ortolog,file="CCArag_ortolog.Rdata")

#Dieta
ortolog_Dieta <- orthologs(genes=Dieta_vs_normal$ID,species="Mus musculus",human=F)
Dieta_vs_normal_ortolog <- rename(Dieta_vs_normal,symbol="ID")
Dieta_vs_normal_ortolog <- merge(ortolog_Dieta,Dieta_vs_normal_ortolog,
                                  by="symbol")
save(Dieta_vs_normal_ortolog,file="Dieta_ortolog.Rdata")

#Hidrodinamica
ortolog_Hidrodinamica <- orthologs(genes=Hidrodinamica_vs_normal$ID,species="Mus musculus",human=F)
Hidrodinamica_vs_normal_ortolog <- rename(Hidrodinamica_vs_normal,symbol="ID")
Hidrodinamica_vs_normal_ortolog <- merge(ortolog_Hidrodinamica,Hidrodinamica_vs_normal_ortolog,
                                 by="symbol")
save(Hidrodinamica_vs_normal_ortolog,file="Hidrodinamica_ortolog.Rdata")

#C.Transformados
ortolog_Transformados <- orthologs(genes=Transformadas_vs_normal$ID,species="Mus musculus",human=F)
Transformados_vs_normal_ortolog <- rename(Transformadas_vs_normal,symbol="ID")
Transformados_vs_normal_ortolog <- merge(ortolog_Transformados,Transformados_vs_normal_ortolog,
                                         by="symbol")
save(Transformados_vs_normal_ortolog,file="Transformadas_ortolog.Rdata")

#Dieta Vs Hidrodinámica
ortolog_DietaVsHidro<- orthologs(genes=Dieta_vs_Hidrodinamica$ID,species="Mus musculus",human=F)
DietaVsHidro_ortolog <- rename(Dieta_vs_Hidrodinamica,symbol="ID")
DietaVsHidro_ortolog <- merge(ortolog_DietaVsHidro,DietaVsHidro_ortolog,
                                         by="symbol")
save(DietaVsHidro_ortolog ,file="DietaVsHidro_ortolog.Rdata")
