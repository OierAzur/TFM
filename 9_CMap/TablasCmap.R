#Trabajo Fin de Máster - Máster en Métodos Computacionales en ciencias
#Autor: Oier Azurmendi Senar
#Tutores_ Silvestre Vicent Cambra y Mikel Hernaez
#Curso 2021-2022

#Objetivo: Preparación de las listas con los nombres de los genes para utilizar como input del CMap.

######################################################################################################
setwd("Documents/TFM/CMap/")

library(babelgene)
load("/Users/oierazur/TFM/DGE/CCARag_vs_normal_filtro1.Rdata")
load("/Users/oierazur/TFM/DGE/Dieta_vs_normal_filtro1.Rdata")
load("/Users/oierazur/TFM/DGE/Hidrodinamica_vs_normalfiltro1.Rdata")
load("/Users/oierazur/TFM/DGE/Transformadas_vs_normalfiltro1.Rdata")
load("/Users/oierazur/TFM/DGE/Dieta_vs_Hidrodinamica.Rdata")

library(babelgene)

#CCARAG vs NORMAL
#UP-regulated
head(CCARag_vs_normal)
CCARag_vs_normal_mvg_up <- CCARag_vs_normal[ order(CCARag_vs_normal$logFC,decreasing=T), ][1:200,1]
ortolog_CCARag_up <- orthologs(genes=CCARag_vs_normal_mvg_up,species="Mus musculus",human=F)
CCARagUP <- ortolog_CCARag_up$human_symbol
write.table(CCARagUP ,file="CCARagUP.txt",row.names = F,col.names = "F",quote = F)
#DOWN-regulated
CCARag_vs_normal_mvg_down <- CCARag_vs_normal[ order(CCARag_vs_normal$logFC), ][1:200,1]
ortolog_CCARag_down <- orthologs(genes=CCARag_vs_normal_mvg_down,species="Mus musculus",human=F)
CCARagDOWN <- ortolog_CCARag_down$human_symbol
write.table(CCARagDOWN ,file="CCARagDOWN.txt",row.names = F,col.names = "F",quote = F)

#Dieta VS Normal
#UP-regulated
Dieta_vs_normal_mvg_up <- Dieta_vs_normal[ order(Dieta_vs_normal$logFC,decreasing=T), ][1:200,1]
ortolog_Dieta_up <- orthologs(genes=Dieta_vs_normal_mvg_up,species="Mus musculus",human=F)
DietaUP <- ortolog_Dieta_up$human_symbol
write.table(DietaUP ,file="DietaUP.txt",row.names = F,col.names = "F",quote = F)
#DOWN-regulated
ortolog_Dieta_down <- orthologs(genes=Dieta_vs_normal_mvg_down,species="Mus musculus",human=F)
Dieta_vs_normal_mvg_down <- Dieta_vs_normal[ order(Dieta_vs_normal$logFC), ][1:200,1]
DietaDOWN <- ortolog_Dieta_down$human_symbol
write.table(DietaDOWN ,file="DietaDOWN.txt",row.names = F,col.names = "F",quote = F)

#Hidrodinamica vs normal
#UP-regulated
head(Hidrodinamica_vs_normal)
Hidrodinamica_vs_normal_mvg_up <- Hidrodinamica_vs_normal[ order(Hidrodinamica_vs_normal$logFC,decreasing=T), ][1:200,1]
ortolog_Hidrodinamica_up <- orthologs(genes=Hidrodinamica_vs_normal_mvg_up,species="Mus musculus",human=F)
HidrodinamicaUP <- ortolog_Hidrodinamica_up$human_symbol
write.table(HidrodinamicaUP ,file="HidrodinamicaUP.txt",row.names = F,col.names = "F",quote = F)
#DOWN-regulated
Hidrodinamica_vs_normal_mvg_down <- Hidrodinamica_vs_normal[ order(Hidrodinamica_vs_normal$logFC), ][1:200,1]
ortolog_Hidrodinamica_down <- orthologs(genes=Hidrodinamica_vs_normal_mvg_down,species="Mus musculus",human=F)
HidrodinamicaDOWN <- ortolog_Hidrodinamica_down$human_symbol
write.table(HidrodinamicaDOWN ,file="HidrodinamicaDOWN.txt",row.names = F,col.names = "F",quote = F)

#Celulas Transformadas vs normal
#UP-regulated
head(Transformadas_vs_normal)
Transformadas_vs_normal_mvg_up <- Transformadas_vs_normal[ order(Transformadas_vs_normal$logFC,decreasing=T), ][1:200,1]
ortolog_Transformadas_up <- orthologs(genes=Transformadas_vs_normal_mvg_up,species="Mus musculus",human=F)
TransformadasUP <- ortolog_Transformadas_up$human_symbol
write.table(TransformadasUP ,file="TransformadasUP.txt",row.names = F,col.names = "F",quote = F)
#DOWN-regulated
Transformadas_vs_normal_mvg_down <-Transformadas_vs_normal[ order(Transformadas_vs_normal$logFC), ][1:200,1]
ortolog_Transformadas_down <- orthologs(genes=Transformadas_vs_normal_mvg_down,species="Mus musculus",human=F)
TransformadasDOWN <- ortolog_Transformadas_down$human_symbol
write.table(TransformadasDOWN ,file="TransformadasDOWN.txt",row.names = F,col.names = "F",quote = F)

#Dieta vs Hidrodinamica
#UP-regulated
head(Dieta_vs_Hidrodinamica)
Dieta_vs_Hidrodinamica_mvg_up <- Dieta_vs_Hidrodinamica[ order(Dieta_vs_Hidrodinamica$logFC,decreasing=T), ][1:200,1]
ortolog_Dieta_vs_Hidrodinamica_up <- orthologs(genes=Dieta_vs_Hidrodinamica_mvg_up,species="Mus musculus",human=F)
DIETAvsHidrodinamicaUP <- ortolog_Dieta_vs_Hidrodinamica_up$human_symbol
write.table(DIETAvsHidrodinamicaUP  ,file="DIETAvsHidrodinamicaUP .txt",row.names = F,col.names = "F",quote = F)
#DOWN-regulated
Dieta_vs_Hidrodinamica_mvg_down <-Dieta_vs_Hidrodinamica[ order(Dieta_vs_Hidrodinamica$logFC), ][1:200,1]
ortolog_Dieta_vs_Hidrodinamica_down<- orthologs(genes=Dieta_vs_Hidrodinamica_mvg_down,species="Mus musculus",human=F)
DIETAvsHidrodinamicaDOWN <- ortolog_Dieta_vs_Hidrodinamica_down$human_symbol
write.table(DIETAvsHidrodinamicaDOWN  ,file="DIETAvsHidrodinamicaDOWN .txt",row.names = F,col.names = "F",quote = F)
