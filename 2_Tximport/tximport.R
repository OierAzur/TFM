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
library(EnsDb.Mmusculus.v79)
txdb <- EnsDb.Mmusculus.v79
k <- keys(txdb,keytype = "GENEID",columns="TXNAME")
df <- AnnotationDbi::select(txdb,keys=k,keytype = "GENEID", columns = "TXNAME")
tx2gene <- df[, 2:1] 
head(tx2gene)

#Tximport object a nivel de gen
txi = tximport(kal_dirs, type = "kallisto", 
               tx2gene=tx2gene,ignoreTxVersion=TRUE)
rownames(s2c)
rownames(txi$counts)
