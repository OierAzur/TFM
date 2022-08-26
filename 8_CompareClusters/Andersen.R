#Trabajo Fin de Máster - Máster en Métodos Computacionales en ciencias
#Autor: Oier Azurmendi Senar
#Tutores_ Silvestre Vicent Cambra y Mikel Hernaez
#Curso 2021-2022

#Objetivo: Observamos la relación de nuestras muestras en los subgrupos publicados por Andersen

######################################################################################################

setwd("Documents/TFM/ClusterComparison/Andersen")

load("/Users/TFM/DGE/ObjetoTxi.Rdata")
load("/Users/TFM/DGE/metadata.Rdata")

#Cargamos los datos de Andersen
library(GEOquery)
library(dplyr)
library(limma)
library(edgeR)
library(AnnotationDbi)
library(RColorBrewer)
library(gplots)
library(EnsDb.Mmusculus.v79)

txdb <- EnsDb.Mmusculus.v79

## change my_id to be the dataset that you want.  Andersen
my_id <- "GSE26566"
gset<- getGEO(my_id,GSEMatrix =TRUE, getGPL=FALSE)
geo_Info <- pData(gset[[1]])

#Cargamos las tablas necesarias
ex<- read.csv2("GSE26566_non-normalized.csv",header=T,sep=";",dec=".")
classifier <- read.csv2("showFullTableHTML.csv",header=T,skip = 1,dec=".")


ex_for_DGE <- ex[,2:170]
ex_DGE <- DGEList(ex_for_DGE ,samples=colnames(ex_for_DGE ))
rownames(ex_DGE$counts) <- ex$Gene.Symbol

#Genes clasificadores
keep <- which(rownames(ex_DGE)%in%classifier$Gene.symbol)
ex_DGE <- ex_DGE[keep,]

#CEACAM - Para clasificas los genes
counts_CEACAM <- cpm(ex_DGE$counts["TMPRSS4",])
sampleInfo=data.frame(vector())
for (i in 1:169){
  if(counts_CEACAM[i,]>3000 ){
    sampleInfo[i,] <- "Poor prognostic" 
  } else{ 
    sampleInfo[i,] <- "Good prognostic"
  }
}

#Consigo mis datos
datos <- txi$counts

rownames(datos) <- mapIds(txdb,
                          keys=rownames(datos),
                          column="SYMBOL",
                          keytype="GENEID",
                          multiVals="first")
library(orthogene)
datos <- orthogene::convert_orthologs(gene_df = datos,
                                          gene_input = "rownames", 
                                          gene_output = "rownames", 
                                          input_species = "mouse",
                                          output_species = "human",
                                          non121_strategy = "drop_both_species") 


unidos <- merge(datos,ex_DGE$counts,by=0)
rownames(unidos) <- unidos$Row.names
unidos <- unidos[,2:185]
nombres <- colnames(ex_DGE$counts)
samples <- c(s2c$sample,nombres)
group <- as.factor(c(s2c$condition,sampleInfo$vector..))

mis_datos <- DGEList(unidos,samples=samples,group=group)

#Normalización
keep <- rowSums(cpm(mis_datos)>2)>= 3
y <- mis_datos[keep, ]
keep <- filterByExpr(y)
y <- y[keep, ]

#counts por millon mayor que 0.5
#Filtrado de los genes poco expresados
treshold <- y$counts>0.5
head(treshold)
#Resumen de los trues en cada muestra
table(rowSums(treshold))
keep <- rowSums(treshold)>=2
summary(keep)

y <- calcNormFactors(y)

#log2 counts per million - Ver la normalización
logcounts <- cpm.DGEList(y,log=T,normalized.lib.sizes = T,prior.count = 5)
boxplot(logcounts,xlab="",ylab="log2 counts per million",las=2,outline=F)
abline(h=median(logcounts),col="blue")

#MDS Multidimensional - Clustering no supervisado
plotMDS(y$counts)

#Variability
var_genes <- apply(logcounts,1,var)
head(var_genes)

#Select 100 more variable genes
select_var <- names(sort(var_genes,decreasing=T))[1:100]
head(select_var)
#Subset logcounts matrix
highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm)
head(highly_variable_lcpm)

#Le aplicamos una matriz de diseño

design <- model.matrix(~ 0 + group)
design

## Make the column names of the design matrix a bit nicer

colnames(design) <- levels(group)
design

#Aplicamos voom

par(mfrow=c(1,1))
v <- voom(y,design,plot = TRUE,normalize.method = "quantile")
names(v)

#Generación del Heatmap
colores <- c("green","blue","purple","orange","yellow","red","grey")
col.cell <- colores[mis_datos$samples$group]
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
heatmap.2(v$E ,col=rev(morecols(100)),trace="none",main="Andersen Cluster Comparison",
          ColSideColors=col.cell,scale="row",na.rm=T,dendrogram = "column")
boxplot(v$E,outline=F)
legend("bottomleft",fill=colores,legend=levels(y$samples$group),cex=0.43)


