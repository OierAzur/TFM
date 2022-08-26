#Trabajo Fin de Máster - Máster en Métodos Computacionales en ciencias
#Autor: Oier Azurmendi Senar
#Tutores_ Silvestre Vicent Cambra y Mikel Hernaez
#Curso 2021-2022

#Objetivo: Análisis de sobrerrepresentación de los KEGG pathways.

######################################################################################################
setwd("Documents/TFM/KEGG/")

load("/Users/oierazur/TFM/GO/CCArag_ortolog.Rdata")
load("/Users/oierazur/TFM/GO/Dieta_ortolog.Rdata.Rdata")
load("/Users/oierazur/TFM/GO/Hidrodinamica_ortolog.Rdata")
load("/Users/oierazur/TFM/GO/TTransformadas_ortolog.Rdata")
load("/Users/oierazur/TFM/GODietaVsHidro_ortolog.Rdata")

library("clusterProfiler")
library("enrichplot")

#CCA (Rag-F1)
genes_CCARag <- as.character(CCARag_vs_normal_ortolog$human_entrez)

#Enrichment pathways
kk_CCARag <- enrichKEGG(gene=genes_CCARag,organism='hsa',pvalueCutoff = 0.05)
dotplot(kk_CCARag,showCategory=15,title="Enriched Pathways CCARag vs Normal")
#Resultados
results_CCARag <- kk_CCARag@result
#Ordener según Ratio de los genes
results_CCARag <- results_CCARag[order(results_CCARag$GeneRatio,decreasing=T),]
head(results_CCARag)
save(results_CCARag,file="pathways_generation_CCARag.Rdata")
write.csv(x=results_CCARag,"pathways_generation_CCARag.csv")

# Dieta
genes_Dieta <- as.character(Dieta_vs_normal_ortolog$human_entrez)
#Enrichment pathways
kk_Dieta <- enrichKEGG(gene=genes_Dieta,organism='hsa',pvalueCutoff = 0.05)
dotplot(kk_Dieta,showCategory=15,title="Enriched Pathways Dieta vs Normal")
#Resultados
results_Dieta<- kk_Dieta@result
#Ordener según Ratio de los genes
results_Dieta <- results_Dieta[order(results_Dieta$GeneRatio,decreasing=T),]
head(results_Dieta)
save(results_Dieta,file="pathways_generation_Dieta.Rdata")
write.csv(x=results_Dieta,"pathways_generation_Dieta.csv")

# Hidrodinamica

genes_Hidrodinamica <- as.character(Hidrodinamica_vs_normal_ortolog$human_entrez)

#Enrichment pathways
kk_Hidrodinamica  <- enrichKEGG(gene=genes_Hidrodinamica ,organism='hsa',pvalueCutoff = 0.05)
dotplot(kk_Hidrodinamica ,showCategory=15,title="Enriched Pathways Hidrodinamica vs Normal")
#Resultados
results_Hidrodinamica  <- kk_Hidrodinamica@result
#Ordener según Ratio de los genes
results_Hidrodinamica  <- results_Hidrodinamica [order(results_Hidrodinamica$GeneRatio,decreasing=T),]
head(results_Hidrodinamica)
save(results_Hidrodinamica,file="pathways_generation_Hidrodinamica.Rdata")
write.csv(x=results_Hidrodinamica,"pathways_generation_Hidrodinamica.csv")

# C. Transformadas
genes_transformados <- as.character(Transformados_vs_normal_ortolog$human_entrez)

#Enrichment pathways
kk_transformados <- enrichKEGG(gene=genes_transformados,organism='hsa',pvalueCutoff = 0.05)
dotplot(kk_transformados,showCategory=15,title="Enriched Pathways Transformadas vs Normal")
#Resultados
results_transformados <- kk_transformados@result
#Ordener según Ratio de los genes
results_transformados <- results_transformados[order(results_transformados$GeneRatio,decreasing=T),]
head(results_transformados)
save(results_transformados,file="pathways_generation_Transformadas.Rdata")
write.csv(x=results_transformados,"pathways_generation_Transformadas.csv")

# Dieta Vs Hidrodinamica
genes_DieHidr <- as.character(Dieta_vs_Hidro_ortolog$human_entrez)

#Enrichment pathways
kk_DieHidr <- enrichKEGG(gene=genes_DieHidr,organism='hsa',pvalueCutoff = 0.05)
dotplot(kk_DieHidr,showCategory=15,title="Enriched Pathways Dieta Vs Hidrodinamica")
results_DieHidr <- kk_DieHidr@result
#Ordener según Ratio de los genes
results_DieHidr <- results_DieHidr[order(results_DieHidr$GeneRatio,decreasing=T),]
head(results_DieHidr)
save(results_DieHidr,file="pathways_generation_DietaVsHidro.Rdata")
write.csv(x=results_DieHidr,"pathways_generation_DietaVsHidro.csv")



