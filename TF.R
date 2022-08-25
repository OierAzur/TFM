BiocManager::install("msigdb")
library(msigdb)
library(ExperimentHub)
library(GSEABase)
eh = ExperimentHub()
query(eh , 'msigdb')
msigdb.hs = getMsigdb(org = 'hs', id = 'SYM', version = '7.4')
msigdb.hs
msigdb.hs = appendKEGG(msigdb.hs)
msigdb.hs
#Acces to the data
length(msigdb.hs)
gs = msigdb.hs[[1000]]
gs
geneIds(gs)
collectionType(gs)
bcCategory(collectionType(gs))
bcSubCategory(collectionType(gs))
description(gs)
details(gs)
#calculate the number of signatures in each category
table(sapply(lapply(msigdb.hs, collectionType), bcCategory))
table(sapply(lapply(msigdb.hs, collectionType), bcSubCategory))
hist(sapply(lapply(msigdb.hs, geneIds), length),
     main = 'MSigDB signature size distribution',
     xlab = 'Signature size')
#Subset collections
listCollections(msigdb.hs)
listSubCollections(msigdb.hs)
subsetCollection(msigdb.hs, 'h')
subsetCollection(msigdb.hs, 'c5')
#Prepare Collections for Limma
library(limma)
Complete <- Complete[,c("Sample1","Sample2","Sample3","Sample4","Sample5",
                        "Sample6","Sample7","Sample8","Sample9","Sample10","Sample11","Sample12",
                        "Sample13","Sample14","Sample15")]
Complete <- apply(Complete,2,as.numeric)
rownames(Complete) <- our_data$human_symbol
#retrieve collections
hallmarks = subsetCollection(msigdb.hs, 'h')
msigdb_ids = geneIds(hallmarks)
#Meto los Ids de los genes que me interesan
fry_indices = ids2indices(msigdb_ids, CCARag_vs_normal_ortolog$human_symbol)
fry_indices[1:2]
head(fry_indices)


