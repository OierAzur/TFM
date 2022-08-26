#Trabajo Fin de Máster - Máster en Métodos Computacionales en ciencias
#Autor: Oier Azurmendi Senar
#Tutores_ Silvestre Vicent Cambra y Mikel Hernaez
#Curso 2021-2022

#Objetivo: Observamos la relación de nuestras muestras en los subgrupos publicados por Sia en 2013.

######################################################################################################

setwd("Documents/TFM/ClusterComparison/Sia")

load("/Users/TFM/DGE/ObjetoTxi.Rdata")
load("/Users/TFM/DGE/metadata.Rdata")
#cargar datos de expresión no normalizados
ex<- read.table("GSE32225_non-normalized.txt",header=T,dec=".",row.names = 1)
classifier <- t(read.table("classifier_sia.txt",sep=" "))

library(illuminaHumanv4.db)
library(GEOquery)
library(EnsDb.Mmusculus.v79)
library(orthogene)
library(RColorBrewer)
library(gplots)

txdb <- EnsDb.Mmusculus.v79

## change my_id to be the dataset that you want. SIA
my_id <- "GSE32225"
gset<- getGEO(my_id,GSEMatrix =TRUE, getGPL=FALSE)

sampleInfo <- pData(gset[[1]])
sampleInfo$`nmf.subclass:ch1`[150:155] <- "control"

# DGE list
y <- DGEList(ex,samples=colnames(ex),group=sampleInfo$`nmf.subclass:ch1`)
rownames(y$counts) <- mapIds(illuminaHumanv4.db,
                        keys=rownames(y$counts),
                        column="SYMBOL",
                        keytype="PROBEID",
                        multiVals="first")
                        
#Cargamos nuestros datos
datos <- txi$counts
rownames(datos) <- mapIds(txdb,
                         keys=rownames(datos),
                         column="SYMBOL",
                         keytype="GENEID",
                         multiVals="first")


mis_datos <- orthogene::convert_orthologs(gene_df = datos,
                                          gene_input = "rownames", 
                                          gene_output = "rownames", 
                                          input_species = "mouse",
                                          output_species = "human",
                                          non121_strategy = "drop_both_species") 
                                          
unidos <- merge(mis_datos,y$counts,by=0)
rownames(unidos) <- unidos$Row.names
unidos <- unidos[,2:171]
nombres <- colnames(y$counts)
samples <- c(s2c$sample,nombres)
group <- as.factor(c("CCA.Rag2_F1", "CCA.Rag2_F1", "CCA.Rag2_F1", "Dieta", "Dieta", "Dieta", 
                     "Hidrodinamica", "Hidrodinamica", "Hidrodinamica","Colangiocitos.transformados",
                     "Colangiocitos.transformados","Colangiocitos.transformados","Colangiocitos.normales","Colangiocitos.normales",
                     "Colangiocitos.normales",sampleInfo$`nmf.subclass:ch1`))

mis_datos <- DGEList(unidos,samples=samples,group=group)

#Normalización
keep <- rowSums(cpm(mis_datos)>2)>= 3
y <- mis_datos[keep, ]
keep <- filterByExpr(y)
y <- y[keep, ]

#counts por millon mayor que 0.5
#Filtrado de los genes poco expresados
miCPM <- cpm(y)
treshold <- miCPM>0.5
head(treshold)

#Resumen de los trues en cada muestra
table(rowSums(treshold))
keep <- rowSums(treshold)>=2
summary(keep)

#Genes clasificadores
classifier <- as.data.frame(classifier)
keep <- which(rownames(y)%in%classifier$V1)
y <- y[keep,]
y$counts[,1:15] <- log(y$counts[,1:15])
y <- calcNormFactors(y)

#log2 counts per million - Ver la normalización
logcounts <- cpm.DGEList(y,log=T,normalized.lib.sizes = T,prior.count = 5)
boxplot(logcounts,xlab="",ylab="log2 counts per million",las=2,outline=F)
abline(h=median(logcounts),col="blue")

#Variability
var_genes <- apply(logcounts,1,var)
head(var_genes)

#Ver los genes más variables
select_var <- names(sort(var_genes,decreasing=T))
head(select_var)
highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm)
head(highly_variable_lcpm)

#MDS Multidimensional - Clustering no supervisado
plotMDS(highly_variable_lcpm)

#Le aplicamos una matriz de diseño

design <- model.matrix(~ 0 + group)
design

## Make the column names of the design matrix a bit nicer

colnames(design) <- levels(group)
design


#Aplicamos voom

par(mfrow=c(1,1))
v <- voom(y,design,plot = TRUE,normalize.method="quantile)
names(v)

#Generación del Heatmap
colores <- c("green","blue","purple","orange","yellow","red","grey","black")
col.cell <- colores[mis_datos$samples$group]
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
heatmap.2(v$E ,col=rev(morecols(100)),trace="none",main="Sia 2013 Cluster Comparison",
          ColSideColors=col.cell,scale="row",na.rm=T,dendrogram = "column")
legend("bottomleft",fill=colores,legend=levels(y$samples$group),cex=0.43)



