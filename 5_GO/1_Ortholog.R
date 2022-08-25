#Trabajo Fin de Máster - Máster en Métodos Computacionales en ciencias
#Autor: Oier Azurmendi Senar
#Tutores_ Silvestre Vicent Cambra y Mikel Hernaez
#Curso 2021-2022

#Objetivo: Obtener los ortógos de humano para poder realizar los análisis de Gene Ontology

######################################################################################################
library(babelgene)

ortolog_CCARag <- orthologs(genes=CCARag_vs_normal$ID,species="Mus musculus",human=F)
CCARag_vs_normal_ortolog <- rename(CCARag_vs_normal,symbol="ID")
CCARag_vs_normal_ortolog <- merge(ortolog_CCARag,CCARag_vs_normal_ortolog,
                                  by="symbol") 

ortolog_Dieta <- orthologs(genes=Dieta_vs_normal$ID,species="Mus musculus",human=F)
Dieta_vs_normal_ortolog <- rename(Dieta_vs_normal,symbol="ID")
Dieta_vs_normal_ortolog <- merge(ortolog_Dieta,Dieta_vs_normal_ortolog,
                                  by="symbol")

ortolog_Hidrodinamica <- orthologs(genes=Hidrodinamica_vs_normal$ID,species="Mus musculus",human=F)
Hidrodinamica_vs_normal_ortolog <- rename(Hidrodinamica_vs_normal,symbol="ID")
Hidrodinamica_vs_normal_ortolog <- merge(ortolog_Hidrodinamica,Hidrodinamica_vs_normal_ortolog,
                                 by="symbol")

ortolog_Transformados <- orthologs(genes=Transformadas_vs_normal$ID,species="Mus musculus",human=F)
Transformados_vs_normal_ortolog <- rename(Transformadas_vs_normal,symbol="ID")
Transformados_vs_normal_ortolog <- merge(ortolog_Transformados,Transformados_vs_normal_ortolog,
                                         by="symbol")
ortolog_DietaVsHidro<- orthologs(genes=Dieta_vs_Hidrodinamica$ID,species="Mus musculus",human=F)
DietaVsHidro_ortolog <- rename(Dieta_vs_Hidrodinamica,symbol="ID")
DietaVsHidro_ortolog <- merge(ortolog_DietaVsHidro,DietaVsHidro_ortolog,
                                         by="symbol")
