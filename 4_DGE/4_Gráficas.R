#Trabajo Fin de Máster - Máster en Métodos Computacionales en ciencias
#Autor: Oier Azurmendi Senar
#Tutores_ Silvestre Vicent Cambra y Mikel Hernaez
#Curso 2021-2022

#Objetivo: Heatmaps y volcano plotsde DGE en cada contraste y expresión de 
#genes seleccionados en todos los contrastes

######################################################################################################
setwd("Users/oierazur/TFM/DGE/")
library(RColorBrewer)
library(gplots)
library(edgeR)
library(dplyr)
load("/Users/oierazur/TFM/DGE/CCARag_vs_normal_filtro1.Rdata")
load("/Users/oierazur/TFM/DGE/Dieta_vs_normal_filtro1.Rdata")
load("/Users/oierazur/TFM/DGE/Hidrodinamica_vs_normalfiltro1.Rdata")
load("/Users/oierazur/TFM/DGE/Transformadas_vs_normalfiltro1.Rdata")
load("/Users/oierazur/TFM/DGE/Dieta_vs_Hidrodinamica.Rdata")
load("/Users/oierazur/TFM/DGE/voom_result.Rdata")

#Heatmaps
lcpm <- cpm(y,log=T)

#CCARag
CCARag_vs_normal.topgenes <- CCARag_vs_normal$ID[1:1000]
i <- which(rownames(y) %in% CCARag_vs_normal.topgenes)
mycol <- colorpanel(1000,"blue","white","red")
lcpm_CCARag <- lcpm[i,1:3]
lcpm_Normal <- lcpm[i,13:15]
lcpmCCAvsNormal <-cbind(lcpm_CCARag,lcpm_Normal)
heatmap.2(lcpmCCAvsNormal, scale="row",
          labRow=rownames(y)[i], labCol=c("CCARag","CCARag","CCARag","Colangiocitos Normales","Colangiocitos Normales","Colangiocitos Normales"), 
          col=mycol, trace="none",density.info="none",
          margin=c(14,7),lhei=c(2,10), dendrogram="column",main="Most variable genes CCARag")

#Dieta
Dieta_vs_normal.topgenes <- Dieta_vs_normal$ID[1:1000]
i <- which(rownames(y) %in% Dieta_vs_normal.topgenes)
lcpm_Dieta <- lcpm[i,4:6]
lcpmDietavsNormal <-cbind(lcpm_Dieta,lcpm_Normal)
heatmap.2(lcpmDietavsNormal, scale="row",
          labRow=rownames(y)[i], labCol = c("Dieta","Dieta","Dieta","Colangiocitos Normales","Colangiocitos Normales","Colangiocitos Normales"),
          col=mycol, trace="none",density.info="none",
          margin=c(14,7),lhei=c(2,10), dendrogram="column",main="Most variable genes Dieta")

#Hidrodinamica
Hidrodinamica_vs_normal.topgenes <- Hidrodinamica_vs_normal$ID[1:999]
i <- which(rownames(y) %in% Hidrodinamica_vs_normal.topgenes)
lcpm_Hidrodinamica <- lcpm[i,7:9]
lcpmHidrodinamicavsNormal <-cbind(lcpm_Hidrodinamica,lcpm_Normal)
heatmap.2(lcpmHidrodinamicavsNormal, scale="row",
          labRow=rownames(y)[i], labCol=c("Hidrodinamica","Hidrodinamica","Hidrodinamica","Colangiocitos Normales","Colangiocitos Normales","Colangiocitos Normales"),
          col=mycol, trace="none",density.info="none",
          margin=c(14,7),lhei=c(2,10), dendrogram="column",main="Most variable genes Hidrodinamica")

#C.Transformados
Transformadas_vs_normal.topgenes <- Transformadas_vs_normal$ID[1:1000]
i <- which(rownames(y) %in% Transformadas_vs_normal.topgenes)
lcpm_Transformadas <- lcpm[i,10:12]
lcpmTransformadasvsNormal <-cbind(lcpm_Transformadas,lcpm_Normal)
heatmap.2(lcpmTransformadasvsNormal, scale="row",
          labRow=rownames(y)[i], labCol = c("Transformadas","Transformadas","Transformadas","Colangiocitos Normales","Colangiocitos Normales","Colangiocitos Normales"),
          col=mycol, trace="none",density.info="none",
          margin=c(14,7),lhei=c(2,10), dendrogram="column",main="Most variable genes Transformadas")

#Dieta vs hidrodinamica
Dieta_vs_Hidrodinamica.topgenes <- Dieta_vs_Hidrodinamica$ID[1:1000]
i <- which(rownames(y) %in% Dieta_vs_Hidrodinamica.topgenes)
lcpm_Dieta <- lcpm[i,4:6]
lcpm_Hidrodinamica <- lcpm[i,7:9]
lcpmDieta_vs_Hidrodinamica <-cbind(lcpm_Dieta,lcpm_Hidrodinamica)
heatmap.2(lcpmDieta_vs_Hidrodinamica, scale="row",
          labRow=rownames(y)[i], labCol = c("Dieta","Dieta","Dieta","Hidrodinamica","Hidrodinamica","Hidrodinamica"),
          col=mycol, trace="none",density.info="none",
          margin=c(14,7),lhei=c(2,10), dendrogram="column",main="Most variable Dieta vs Hidrodinamica")


#Volcanoplot
library(EnhancedVolcano)
volcanoplot(fit.cont,coef = 1,highlight=11,names=rownames(fit.cont$coefficients),main="CCARag Volcano Plot")
Rag_p<- EnhancedVolcano(CCARag_vs_normal,
                        title = 'CCA(Rag-F1) vs Normal',
                lab =CCARag_vs_normal[,1],
                x = 'logFC',
                y = 'P.Value')
Rag_p
volcanoplot(fit.cont,coef = 2,highlight=11,names=rownames(fit.cont$coefficients),main="Dieta Volcano Plot")
Dieta_p<- EnhancedVolcano(Dieta_vs_normal,
                title = 'Dieta vs Normal',
                lab =Dieta_vs_normal[,1],
                x = 'logFC',
                y = 'P.Value')
Dieta_p
volcanoplot(fit.cont,coef = 3,highlight=11,names=rownames(fit.cont$coefficients),main="Hidrodinamica Volcano Plot")
Hidro_p <- EnhancedVolcano(Hidrodinamica_vs_normal,
                title = 'Hidrodinamica vs Normal',
                lab =Hidrodinamica_vs_normal[,1],
                x = 'logFC',
                y = 'P.Value')
Hidro_p 
 volcanoplot(fit.cont,coef = 4,highlight=11,names=rownames(fit.cont$coefficients),main="Transformadas Volcano Plot")
Trans_p <- EnhancedVolcano(Transformadas_vs_normal,
                title = 'Transformados vs Normal',
                lab =Transformadas_vs_normal[,1],
                x = 'logFC',
                y = 'P.Value')
Trans_p
volcanoplot(fit.cont,coef = 5,highlight=11,names=rownames(fit.cont$coefficients),main="Transformadas Volcano Plot")
DH_p <- EnhancedVolcano(Dieta_vs_Hidrodinamica,
                title = 'Dieta Vs Hidrodinamica',
                lab =Dieta_vs_Hidrodinamica[,1],
                x = 'logFC',
                y = 'P.Value')
DH_p 

#Gráficas despues de testar la expresión diferencial
#Despues para cada uno de los contrastes se puede hacer esto:
#Mirar valores de muestras individuales

par(mfrow=c(2,2))
#Stripchart for each DGE gene
Alb <- stripchart(v$E["Alb",]~v$targets$group,vertical=T,las=2,cex.axis=0.52,pch=16,col=1:5,method="jitter",
           ylab="Normalised log2 expression",main="Alb")
           
Trf <- stripchart(v$E["Trf",]~v$targets$group,vertical=T,las=2,cex.axis=0.52,pch=16,col=1:5,method="jitter",
           ylab="Normalised log2 expression",main="Trf")
Rps3a3 <- stripchart(v$E["Rps3a3",]~v$targets$group,vertical=T,las=2,cex.axis=0.52,pch=16,col=1:5,method="jitter",
                     ylab="Normalised log2 expression",main="Rps3a3")
Slc25a48 <- stripchart(v$E["Slc25a48",]~v$targets$group,vertical=T,las=2,cex.axis=0.52,pch=16,col=1:5,method="jitter",
           ylab="Normalised log2 expression",main="Slc25a48")
           
#Para el gen que se quiera:
#gene_name<- "nombre gen"
#stripchart(v$E[gene_name,]~v$targets$group,vertical=T,las=2,cex.axis=0.52,pch=16,col=1:5,method="jitter",
           ylab="Normalised log2 expression",main=gene_name)           
          
