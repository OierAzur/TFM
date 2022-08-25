#Trabajo Fin de Máster - Máster en Métodos Computacionales en ciencias
#Autor: Oier Azurmendi Senar
#Tutores_ Silvestre Vicent Cambra y Mikel Hernaez
#Curso 2021-2022

#Objetivo: Importar la abundancia, la cuantificación estimada y la longitud a nivel de tránscrito para
#poder realizar los análisis a nivel de gen.

######################################################################################################
library(AnnotationDbi)
library(EnsDb.Mmusculus.v79)
setwd("Documents/TFM/Kallisto_complete_Results/")

base_dir <- "."
sample_ids <- dir(base_dir,"_")
sample_ids
kal_dirs <- file.path( sample_ids, "abundance.h5")
kal_dirs
s2c = read.table(file.path(base_dir,"design.txt"), header=TRUE,
                 stringsAsFactors=FALSE,sep ="-")
s2c
s2c$path=kal_dirs
print(s2c)

#tximportData Object

txdb <- EnsDb.Mmusculus.v79
k <- keys(txdb,keytype = "GENEID",columns="TXNAME")
df <- AnnotationDbi::select(txdb,keys=k,keytype = "GENEID", columns = "TXNAME")
tx2gene <- df[, 2:1] 
head(tx2gene)

#Tximport object a nivel de gen
txi = tximport(kal_dirs, type = "kallisto", 
               tx2gene=tx2gene,ignoreTxVersion=TRUE)

#Comprobación
rownames(s2c)
rownames(txi$counts)

save(txi,file="ObjetoTxi.Rdata")
