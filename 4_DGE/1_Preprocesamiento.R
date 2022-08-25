#Trabajo Fin de Máster - Máster en Métodos Computacionales en ciencias
#Autor: Oier Azurmendi Senar
#Tutores_ Silvestre Vicent Cambra y Mikel Hernaez
#Curso 2021-2022

#Objetivo: Realizar el preprocesamiento de datos para realizar el análisis de expresión diferencial 
#de los genes

######################################################################################################
setwd("Documents/TFM/DGE/")
library(limma)
library(edgeR) 
library(AnnotationDbi)

load("/Users/oierazur/TFM/DGE/ObjetoTxi.Rdata")

#Crear objeto DGE
y <- DGEList(txi$counts,samples=s2c$sample,group=s2c$condition)

#Comprobación
colnames(y$counts)

#Eliminar los genes poco expresados
keep <- filterByExpr(y)
y <- y[keep, ]

#Observar los resultados de cuantificación y filtrar
miCPM <- cpm(y)
lcpm <- cpm(y,log=T)
#counts por millon mayor que 0.5
treshold <- miCPM>0.5
head(treshold)
#Resumen de los trues en cada muestra
table(rowSums(treshold))
keep <- rowSums(treshold)>=2
summary(keep)
y <- y[keep,keep.lib.sizes=F]

#Normalizacion de las distribuciones de expresión
y <- calcNormFactors(y)

#Graficas
plot(miCPM[,1],y$counts[,1])

#Control de calidad
y$samples$lib.size
barplot(y$samples$lib.size/1e06,names=colnames(y),las=2,ann=F,cex.names=0.75)
mtext(side=1,text="Muestras",line=4)
mtext(side=2,text="Library size (millions)",line=3)
title("Barplot de los tamaños de las librerias")

#log2 counts por million - ver si están normalizados
logcounts <- cpm(y,log=T)
boxplot(logcounts,xlab="",ylab="log2 counts per million",las=2)
abline(h=median(logcounts),col="blue")

#Clustering no supervisado de las muestras MDS Multidimensional
#Compruebo si los grupos estan bien
levels(y$samples$group)
col.cell <- c("green","blue","purple","orange","yellow")[y$samples$group]
data.frame(y$samples$group,col.cell)
plotMDS(y,col=col.cell,pch = 1)
legend("bottomright",fill=c("green","blue","purple","orange","yellow"),legend=levels(y$samples$group),cex=0.35)
title("Por grupos")

#Hierarchical Clustering con heatmaps
var_genes <- apply(logcounts,1,var)
head(var_genes)

#Select 100 more variable genes
select_var <- names(sort(var_genes,decreasing=T))[1:100]
head(select_var)

#Subset logcounts matrix
highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm)
head(highly_variable_lcpm)

#Heatmap
library(RColorBrewer)
library(gplots)
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
col.cell <- c("green","blue","purple","orange","yellow")[y$samples$group]
rownames(highly_variable_lcpm) <- mapIds(txdb,
                                         keys=rownames(highly_variable_lcpm),
                                         column="SYMBOL",
                                         keytype="GENEID",
                                         multiVals="first")
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none",main="Top 50 most variable genes",
          ColSideColors=col.cell,scale="row")
