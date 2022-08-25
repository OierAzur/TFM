#Trabajo Fin de Máster - Máster en Métodos Computacionales en ciencias
#Autor: Oier Azurmendi Senar
#Tutores_ Silvestre Vicent Cambra y Mikel Hernaez
#Curso 2021-2022

#Objetivo: Obtención de las tablas de DGE para cada condición.

######################################################################################################
load("/Users/oierazur/TFM/DGE/fit_DGE.Rdata")
library(limma)
library(tidyverse)
library(dplyr)

#Obtención de las tablas:

#Filtros posibles
#filtro0 <-B>0
#filtr1<-B>0&(logFC>0.75|logFC< -0.75)
#filtro2<-B>0&(logFC>1|logFC< -1)

#CCA (Rag-F1):
CCARag_vs_normal<- topTable(fit.cont,coef=1,number=Inf,sort.by = "B")
CCARag_vs_normal<-dplyr::filter(CCARag_vs_normal, B >0&(logFC>1|logFC< -1))
CCARag_vs_normal_up <- CCARag_vs_normal[CCARag_vs_normal$logFC>0,]
CCARag_vs_normal_down <- CCARag_vs_normal[CCARag_vs_normal$logFC<0,]
CCARag_vs_normal_most_up <- CCARag_vs_normal[ order(CCARag_vs_normal$logFC,decreasing=T), ][1:200,0]
CCARag_vs_normal_most_down <- CCARag_vs_normal[ order(CCARag_vs_normal$logFC), ][1:200,1]
save(CCARag_vs_normal,file="CCARag_vs_normal_filtro1.Rdata")
write.csv(x=CCARag_vs_normal,"CCARag_vs_normal_filtro1.csv")

#Dieta(Rag-F1):
Dieta_vs_normal<- topTable(fit.cont,coef=2,number=Inf,sort.by = "B")
Dieta_vs_normal<-dplyr::filter(Dieta_vs_normal, B >0&(logFC>1|logFC< -1)) 
Dieta_vs_normal_up<- Dieta_vs_normal[Dieta_vs_normal$logFC>0,]
Dieta_vs_normal_down <- Dieta_vs_normal[Dieta_vs_normal$logFC<0,]
Dieta_vs_normal_most_up <- Dieta_vs_normal[ order(Dieta_vs_normal$logFC,decreasing=T), ][1:200,0]
Dieta_vs_normall_most_down <- Dieta_vs_normal[ order(Dieta_vs_normal$logFC), ][1:200,1]
save(Dieta_vs_normal,file="Dieta_vs_normal_filtro1.Rdata")
write.csv(x=Dieta_vs_normal,"Dieta_vs_normal_filtro1.csv")

#Hidrodinamica
Hidrodinamica_vs_normal<- topTable(fit.cont,coef=3,number=Inf,sort.by = "B")
Hidrodinamica_vs_normal<-dplyr::filter(Hidrodinamica_vs_normal, B >0&(logFC>1|logFC< -1)) 
Hidrodinamica_vs_normal_up <- Hidrodinamica_vs_normal[Hidrodinamica_vs_normal$logFC>0,]
Hidrodinamica_vs_normal_down <- Hidrodinamica_vs_normal[Hidrodinamica_vs_normal$logFC<0,]
Hidrodinamica_vs_normal_most_up <- Hidrodinamica_vs_normal[ order(Hidrodinamica_vs_normal$logFC,decreasing=T), ][1:200,0]
Hidrodinamica_vs_normal_most_down <- Hidrodinamica_vs_normal[ order(Hidrodinamica_vs_normal$logFC), ][1:200,1]
save(Hidrodinamica_vs_normal,file="Hidrodinamica_vs_normalfiltro1.Rdata")
write.csv(x=Hidrodinamica_vs_normal,"Hidrodinamica_vs_normalfiltro1.csv")

#Colangiocitos Transformados
Transformadas_vs_normal<- topTable(fit.cont,coef=4,number=Inf,sort.by = "B")
Transformadas_vs_normal<-dplyr::filter(Transformadas_vs_normal, B >0&(logFC>1|logFC< -1)) 
Transformadas_vs_normal_up <- Transformadas_vs_normal[Transformadas_vs_normal$logFC>0,]
Transformadas_vs_normal_down <- Transformadas_vs_normal[Transformadas_vs_normal$logFC<0,]
Transformadas_vs_normal_most_up <- Transformadas_vs_normal[ order(Transformadas_vs_normal$logFC,decreasing=T), ][1:200,0]
Transformadas_vs_normal_most_down <- Transformadas_vs_normal[ order(Transformadas_vs_normal$logFC), ][1:200,1]
save(Transformadas_vs_normal,file="Transformadas_vs_normalfiltro1.Rdata")
write.csv(x=Transformadas_vs_normal,"Transformadas_vs_normalfiltro1.csv")

#Dieta Vs Hidrodinamica
Dieta_vs_Hidrodinamica<- topTable(fit.cont,coef=5,number=Inf,sort.by = "B")
Dieta_vs_Hidrodinamica<-dplyr::filter(Dieta_vs_Hidrodinamica, B >0&(logFC>1|logFC< -1)) 
Dieta_vs_Hidrodinamica_up <- Dieta_vs_Hidrodinamica[Dieta_vs_Hidrodinamica$logFC>0,]
Dieta_vs_Hidrodinamica_down <- Dieta_vs_Hidrodinamica[Dieta_vs_Hidrodinamica$logFC<0,]
Dieta_vs_Hidrodinamica_most_up <- Dieta_vs_Hidrodinamica[ order(TDieta_vs_Hidrodinamica$logFC,decreasing=T), ][1:200,0]
Dieta_vs_Hidrodinamica_most_down <-Dieta_vs_Hidrodinamica[ order(Dieta_vs_Hidrodinamica$logFC), ][1:200,1]
save(Dieta_vs_Hidrodinamica,file="Dieta_vs_Hidrodinamica.Rdata")
write.csv(x=Dieta_vs_Hidrodinamica,"Dieta_vs_Hidrodinamica.csv")

#Tabla General 
tabla_completa<- topTable(fit.cont,coef=NULL,number=Inf,sort.by = "B",resort.by="logFC")
save(tabla_completa,file="tabla_completa_DEG.Rdata")
write.csv(x=tabla_completa,"tabla_completa_DEG.csv")
