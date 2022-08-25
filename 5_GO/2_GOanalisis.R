#Trabajo Fin de Máster - Máster en Métodos Computacionales en ciencias
#Autor: Oier Azurmendi Senar
#Tutores_ Silvestre Vicent Cambra y Mikel Hernaez
#Curso 2021-2022

#Objetivo: Análisis de Gene Ontology y sobrerrepresentación de Gene Ontology.

######################################################################################################
setwd("Documents/TFM/GO/")

load("/Users/oierazur/TFM/GO/CCArag_ortolog.Rdata")
load("/Users/oierazur/TFM/GO/Dieta_ortolog.Rdata.Rdata")
load("/Users/oierazur/TFM/GO/Hidrodinamica_ortolog.Rdata")
load("/Users/oierazur/TFM/GO/TTransformadas_ortolog.Rdata")
load("/Users/oierazur/TFM/GODietaVsHidro_ortolog.Rdata")

library("clusterProfiler")
library("enrichplot")
library(org.Hs.eg.db)
OrgDB <- org.Hs.eg.db

#CCA (Rag-F1)

genes <- as.character(CCARag_vs_normal_ortolog$human_entrez)
ontologia<- "CC" #También se han analizado #"MF" y "BP".

ggo <- clusterProfiler::groupGO(gene=genes,
                                OrgDb=OrgDB ,ont=ontología,level=3,readable=T)
                                
#Observación de los datos
head(as.data.frame(ggo)[,-5])

#Obtención Tablas
save(ggo,file="GO_CC_CCARag.Rdata") 
#También: save(ggo,file="GO_MF_CCARag.Rdata") y save(ggo,file="GO_BP_CCARag.Rdata")
write.csv(x=ggo,"GO_CC_CCARag.csv") 
#También: write.csv(x=ggo,"GO_MF_CCARag.csv") y write.csv(x=ggo,"GO_BP_CCARag.csv")

#Gráficas
barplot(ggo,drop=T,showCategory=20,vertex.label.cex=0.8)
upsetplot(ggo)

#Test de sobre representación GO
ontologia<- "ALL" #También se han analizado #"CC, #"MF" y "BP".
ego <- clusterProfiler::enrichGO(gene=genes,OrgDb=OrgDB,ont = ontologia,
                                 pAdjustMethod="BH",pvalueCutoff = 0.01,qvalueCutoff = 0.01,
                                 readable=T)
upsetplot(ego)
dotplot(ego,showCategory=15,title="Enriched ORA CCARag vs Normal")
barplot(ego,showCategory=15)
save(ego,file="GO_OVER_ALL_CCARag.Rdata") 
#También: save(ego,file="GO_OVER_CC_CCARag.Rdata"),save(ego,file="GO_OVER_MF_CCARag.Rdata") y save(ego,file="GO_OVER_BP_CCARag.Rdata")
write.csv(x=ego,"GO_OVER_ALL_CCARag.csv") 
#También: write.csv(x=ego,"GO_OVER_CC_CCARag.csv"),write.csv(x=ego,"GO_OVER_MF_CCARag.csv") y write.csv(x=ego,"GO_OVER_BP_CCARag.csv")

#Obtención de la lista de genes
geneList <- as.list(CCARag_vs_normal_ortolog$logFC)
names(geneList) <- as.character(unique(CCARag_vs_normal_ortolog$human_entrez))
geneList <- sort(unlist(geneList),decreasing=T)
head(geneList)

#GSEA ANALYSIS
gsea_CCARag<- gseGO(geneList,ont=ontologia,OrgDb = OrgDB, minGSSize = 15,
                maxGSSize = 500,eps = 0, seed=T,pvalueCutoff = 0.25)
gseaplot(gsea_CCARag,geneSetID = 2)

#Dieta
genes <- as.character(Dieta_vs_normal_ortolog$human_entrez)
ontologia<- "CC" #También se han analizado #"MF" y "BP".

ggo <- clusterProfiler::groupGO(gene=genes,
                                OrgDb=OrgDB ,ont=ontología,level=3,readable=T)
                                
#Observación de los datos
head(as.data.frame(ggo)[,-5])

#Obtención Tablas
save(ggo,file="GO_CC_Dieta.Rdata") 
#También: save(ggo,file="GO_MF_Dieta.Rdata") y save(ggo,file="GO_BP_Dieta.Rdata")
write.csv(x=ggo,"GO_CC_Dieta.csv") 
#También: write.csv(x=ggo,"GO_MF_Dieta.csv") y write.csv(x=ggo,"GO_BP_Dieta.csv")

#Gráficas
barplot(ggo,drop=T,showCategory=20,vertex.label.cex=0.8)
upsetplot(ggo)

#Test de sobre representación GO
ontologia<- "ALL" #También se han analizado #"CC, #"MF" y "BP".
ego <- clusterProfiler::enrichGO(gene=genes,OrgDb=OrgDB,ont = ontologia,
                                 pAdjustMethod="BH",pvalueCutoff = 0.01,qvalueCutoff = 0.01,
                                 readable=T)
upsetplot(ego)
dotplot(ego,showCategory=15,title="Enriched ORA Dieta vs Normal")
barplot(ego,showCategory=15)
save(ego,file="GO_OVER_ALL_Dieta.Rdata") 
#También: save(ego,file="GO_OVER_CC_Dieta.Rdata"),save(ego,file="GO_OVER_MF_Dieta.Rdata") y save(ego,file="GO_OVER_BP_Dieta.Rdata")
write.csv(x=ego,"GO_OVER_ALL_CCARag.csv") 
#También: write.csv(x=ego,"GO_OVER_CC_Dieta.csv"),write.csv(x=ego,"GO_OVER_MF_Dieta.csv") y write.csv(x=ego,"GO_OVER_BP_Dieta.csv")

#Obtención de la lista de genes
geneList <- as.list(Dieta_vs_normal_ortolog$logFC)
names(geneList) <- as.character(unique(Dieta_vs_normal_ortolog$human_entrez))
geneList <- sort(unlist(geneList),decreasing=T)
head(geneList)

#GSEA ANALYSIS
gsea_Dieta<- gseGO(geneList,ont=ontologia,OrgDb = OrgDB, minGSSize = 15,
                maxGSSize = 500,eps = 0, seed=T,pvalueCutoff = 0.25)
gseaplot(gsea_Dieta,geneSetID = 2)

#Hidrodinamica
genes <-  as.character(Hidrodinamica_vs_normal_ortolog$human_entrez)
ontologia<- "CC" #También se han analizado #"MF" y "BP".

ggo <- clusterProfiler::groupGO(gene=genes,
                                OrgDb=OrgDB ,ont=ontología,level=3,readable=T)
                                
#Observación de los datos
head(as.data.frame(ggo)[,-5])

#Obtención Tablas
save(ggo,file="GO_CC_Hidrodinamica.Rdata") 
#También: save(ggo,file="GO_MF_Hidrodinamica.Rdata") y save(ggo,file="GO_BP_Hidrodinamica.Rdata")
write.csv(x=ggo,"GO_CC_Dieta.csv") 
#También: write.csv(x=ggo,"GO_MF_Hidrodinamica.csv") y write.csv(x=ggo,"GO_BP_Hidrodinamica.csv")

#Gráficas
barplot(ggo,drop=T,showCategory=20,vertex.label.cex=0.8)
upsetplot(ggo)

#Test de sobre representación GO
ontologia<- "ALL" #También se han analizado #"CC, #"MF" y "BP".
ego <- clusterProfiler::enrichGO(gene=genes,OrgDb=OrgDB,ont = ontologia,
                                 pAdjustMethod="BH",pvalueCutoff = 0.01,qvalueCutoff = 0.01,
                                 readable=T)
upsetplot(ego)
dotplot(ego,showCategory=15,title="Enriched ORA Hidrodinamica vs Normal")
barplot(ego,showCategory=15)
save(ego,file="GO_OVER_ALL_Hidrodinamica.Rdata") 
#También: save(ego,file="GO_OVER_CC_Hidrodinamica.Rdata"),save(ego,file="GO_OVER_MF_Hidrodinamica.Rdata") y save(ego,file="GO_OVER_BP_Hidrodinamica.Rdata")
write.csv(x=ego,"GO_OVER_ALL_CCARag.csv") 
#También: write.csv(x=ego,"GO_OVER_CC_Hidrodinamica.csv"),write.csv(x=ego,"GO_OVER_MF_Hidrodinamica.csv") y write.csv(x=ego,"GO_OVER_BP_Hidrodinamica.csv")

#Obtención de la lista de genes
geneList <- as.list(Hidrodinamica_vs_normal_ortolog$logFC)
names(geneList) <- as.character(unique(Hidrodinamica_vs_normal_ortolog$human_entrez))
geneList <- sort(unlist(geneList),decreasing=T)
head(geneList)

#GSEA ANALYSIS
gsea_Hidrodinamica<- gseGO(geneList,ont=ontologia,OrgDb = OrgDB, minGSSize = 15,
                maxGSSize = 500,eps = 0, seed=T,pvalueCutoff = 0.25)
gseaplot(gsea_Hidrodinamica,geneSetID = 2)

#C. Transformados
genes <-  as.character(Transformados_vs_normal_ortolog$human_entrez)
ontologia<- "CC" #También se han analizado #"MF" y "BP".

ggo <- clusterProfiler::groupGO(gene=genes,
                                OrgDb=OrgDB ,ont=ontología,level=3,readable=T)
                                
#Observación de los datos
head(as.data.frame(ggo)[,-5])

#Obtención Tablas
save(ggo,file="GO_CC_Transformados.Rdata") 
#También: save(ggo,file="GO_MF_Transformados.Rdata") y save(ggo,file="GO_BP_Transformados.Rdata")
write.csv(x=ggo,"GO_CC_Dieta.csv") 
#También: write.csv(x=ggo,"GO_MF_Transformados.csv") y write.csv(x=ggo,"GO_BP_Transformados.csv")

#Gráficas
barplot(ggo,drop=T,showCategory=20,vertex.label.cex=0.8)
upsetplot(ggo)

#Test de sobre representación GO
ontologia<- "ALL" #También se han analizado #"CC, #"MF" y "BP".
ego <- clusterProfiler::enrichGO(gene=genes,OrgDb=OrgDB,ont = ontologia,
                                 pAdjustMethod="BH",pvalueCutoff = 0.01,qvalueCutoff = 0.01,
                                 readable=T)
upsetplot(ego)
dotplot(ego,showCategory=15,title="Enriched ORA C.Transformados vs Normal")
barplot(ego,showCategory=15)
save(ego,file="GO_OVER_ALL_Transformados.Rdata") 
#También: save(ego,file="GO_OVER_CC_Transformados.Rdata"),save(ego,file="GO_OVER_MF_Transformados.Rdata") y save(ego,file="GO_OVER_BP_Transformados.Rdata")
write.csv(x=ego,"GO_OVER_ALL_CCARag.csv") 
#También: write.csv(x=ego,"GO_OVER_CC_Transformados.csv"),write.csv(x=ego,"GO_OVER_MF_Transformados.csv") y write.csv(x=ego,"GO_OVER_BP_Transformados.csv")

#Obtención de la lista de genes
geneList <- as.list(Transformados_vs_normal_ortolog$logFC)
names(geneList) <- as.character(unique(Transformados_vs_normal_ortolog$human_entrez))
geneList <- sort(unlist(geneList),decreasing=T)
head(geneList)

#GSEA ANALYSIS
gsea_Transformados<- gseGO(geneList,ont=ontologia,OrgDb = OrgDB, minGSSize = 15,
                maxGSSize = 500,eps = 0, seed=T,pvalueCutoff = 0.25)
gseaplot(gsea_Transformados,geneSetID = 2)

#Dieta vs Hidrodinamica
genes <-  as.character(DietaVsHidro_ortolog$human_entrez)
ontologia<- "CC" #También se han analizado #"MF" y "BP".

ggo <- clusterProfiler::groupGO(gene=genes,
                                OrgDb=OrgDB ,ont=ontología,level=3,readable=T)
                                
#Observación de los datos
head(as.data.frame(ggo)[,-5])

#Obtención Tablas
save(ggo,file="GO_CC_DietaVsHidro.Rdata") 
#También: save(ggo,file="GO_MF_DietaVsHidro.Rdata") y save(ggo,file="GO_BP_DietaVsHidro.Rdata")
write.csv(x=ggo,"GO_CC_Dieta.csv") 
#También: write.csv(x=ggo,"GO_MF_DietaVsHidro.csv") y write.csv(x=ggo,"GO_BP_DietaVsHidro.csv")

#Gráficas
barplot(ggo,drop=T,showCategory=20,vertex.label.cex=0.8)
upsetplot(ggo)

#Test de sobre representación GO
ontologia<- "ALL" #También se han analizado #"CC, #"MF" y "BP".
ego <- clusterProfiler::enrichGO(gene=genes,OrgDb=OrgDB,ont = ontologia,
                                 pAdjustMethod="BH",pvalueCutoff = 0.01,qvalueCutoff = 0.01,
                                 readable=T)
upsetplot(ego)
dotplot(ego,showCategory=15,title="Enriched ORA Dieta Vs Hidrodinamica")
barplot(ego,showCategory=15)
save(ego,file="GO_OVER_ALL_DietaVsHidro.Rdata") 
#También: save(ego,file="GO_OVER_CC_DietaVsHidro.Rdata"),save(ego,file="GO_OVER_MF_DietaVsHidro.Rdata") y save(ego,file="GO_OVER_BP_DietaVsHidro.Rdata")
write.csv(x=ego,"GO_OVER_ALL_CCARag.csv") 
#También: write.csv(x=ego,"GO_OVER_CC_DietaVsHidro.csv"),write.csv(x=ego,"GO_OVER_MF_DietaVsHidro.csv") y write.csv(x=ego,"GO_OVER_BP_DietaVsHidro.csv")

#Obtención de la lista de genes
geneList <- as.list(DietaVsHidro_ortolog$logFC)
names(geneList) <- as.character(unique(DietaVsHidro_ortolog$human_entrez))
geneList <- sort(unlist(geneList),decreasing=T)
head(geneList)

#GSEA ANALYSIS
gsea_DietaVsHidro<- gseGO(geneList,ont=ontologia,OrgDb = OrgDB, minGSSize = 15,
                maxGSSize = 500,eps = 0, seed=T,pvalueCutoff = 0.25)
gseaplot(gsea_DietaVsHidro,geneSetID = 2)
