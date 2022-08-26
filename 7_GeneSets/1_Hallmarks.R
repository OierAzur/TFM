#Trabajo Fin de Máster - Máster en Métodos Computacionales en ciencias
#Autor: Oier Azurmendi Senar
#Tutores_ Silvestre Vicent Cambra y Mikel Hernaez
#Curso 2021-2022

#Objetivo: Análisis de los Hallmarks de los conjuntos de genes en cada condición

######################################################################################################
setwd("Documents/TFM/GeneSet/Hallmarks")

load("/Users/oierazur/TFM/DGE/voom_result.Rdata")
load("/Users/oierazur/TFM/DGE/contraste.Rdata".Rdata")
load("/Users/oierazur/TFM/GO/CCArag_ortolog.Rdata")
load("/Users/oierazur/TFM/GO/Dieta_ortolog.Rdata.Rdata")
load("/Users/oierazur/TFM/GO/Hidrodinamica_ortolog.Rdata")
load("/Users/oierazur/TFM/GO/TTransformadas_ortolog.Rdata")
load("/Users/oierazur/TFM/GO/DietaVsHidro_ortolog.Rdata")

#MsigDb data
load("/Users/oierazur/TFM/GeneSet/human_H_v5p2.rdata")


#CCARag

idx <- ids2indices(Hs.H,id=CCARag_vs_normal_ortolog$human_entrez)
cam.CCARagsNormal <- camera(v,idx,design,contrast=cont.matrix[,1])
table(cam.CCARagsNormal$FDR < 0.05)
#Ordenado
cam.CCARagsNormal[order(cam.CCARagsNormal$PValue),]
head(cam.CCARagsNormal,5)
#Up-regulados
cam.CCARagsNormal_UP <- subset(x=cam.CCARagsNormal,subset=Direction=="Up")
cam.CCARagsNormal_UP[order(cam.CCARagsNormal_UP$PValue),]
head(cam.CCARagsNormal_UP,5)
#Down-Regulados
cam.CCARagsNormal_Down <- subset(x=cam.CCARagsNormal,subset=Direction=="Down")
cam.CCARagsNormal_Down[order(cam.CCARagsNormal_Down$PValue),]
head(cam.CCARagsNormal_Down,5)
            
save(cam.CCARagsNormal,file="Hallmark_GENESET_CCARag.Rdata")
write.csv(x=cam.CCARagsNormal,"Hallmark_GENESET_CCARag.csv")
save(cam.CCARagsNormal,file="Hallmark_GENESET_CCARag_UP.Rdata")
write.csv(x=cam.CCARagsNormal,"Hallmark_GENESET_CCARag_UP.csv")
save(cam.CCARagsNormal,file="Hallmark_GENESET_CCARag_Down.Rdata")
write.csv(x=cam.CCARagsNormal,"Hallmark_GENESET_CCARag_Down.csv")

#Dieta

idx <- ids2indices(Hs.H,id=Dieta_vs_normal_ortolog$human_entrez)
cam.Dieta_vs_normal <- camera(v,idx,design,contrast=cont.matrix[,2])
table(cam.Dieta_vs_normal$FDR < 0.05)
#Ordenado
cam.Dieta_vs_normal[order(cam.Dieta_vs_normal$PValue),]
head(cam.Dieta_vs_normal,5)
#Up-regulados
cam.Dieta_vs_normal_UP <- subset(x=cam.Dieta_vs_normal,subset=Direction=="Up")
cam.Dieta_vs_normal_UP[order(cam.Dieta_vs_normal_UP$PValue),]
head(cam.Dieta_vs_normal_UP,6)
#Down-regulados
cam.Dieta_vs_normal_Down <- subset(x=cam.Dieta_vs_normal,subset=Direction=="Down")
cam.Dieta_vs_normal_Down[order(cam.Dieta_vs_normal_Down$PValue),]
head(cam.Dieta_vs_normal_Down,6)
            
save(cam.Dieta_vs_normal,file="Hallmark_GENESET_DIETA.Rdata")
write.csv(x=cam.Dieta_vs_normal,"Hallmark_GENESET_DIETA.csv")
save(cam.Dieta_vs_normal_UP,file="Hallmark_GENESET_DIETA_UP.Rdata")
write.csv(x=cam.Dieta_vs_normal_UP,"Hallmark_GENESET_DIETA_UP.csv")
save(cam.Dieta_vs_normal_Down,file="Hallmark_GENESET_DIETA_Down.Rdata")
write.csv(x=cam.Dieta_vs_normal_Down,"Hallmark_GENESET_DIETA_Down.csv")


#Hidrodinamica

idx <- ids2indices(Hs.H,id=Hidrodinamica_vs_normal_ortolog$human_entrez)
cam.Hidrodinamica_vs_normal <- camera(v,idx,design,contrast=cont.matrix[,3])
#Ordenado
cam.Hidrodinamica_vs_normal[order(cam.Hidrodinamica_vs_normal$PValue),]
head(cam.Hidrodinamica_vs_normal,5)
#Up-regulados
cam.Hidrodinamica_vs_normal_UP <- subset(x=cam.Hidrodinamica_vs_normal,subset=Direction=="Up")
cam.Hidrodinamica_vs_normal_UP[order(cam.Hidrodinamica_vs_normal_UP$PValue),]
head(cam.Hidrodinamica_vs_normal_UP,5)
#Down.regulados
cam.Hidrodinamica_vs_normal_Down <- subset(x=cam.Hidrodinamica_vs_normal,subset=Direction=="Down")
cam.Hidrodinamica_vs_normal_Down[order(cam.Hidrodinamica_vs_normal_Down$PValue),]
head(cam.Hidrodinamica_vs_normal_Down,5)
            
save(cam.Hidrodinamica_vs_normal,file="Hallmark_GENESET_HIDRODINAMICA.Rdata")
write.csv(x=cam.Dieta_vs_normal,"Hallmark_GENESET_HIDRODINAMICA.csv")
save(cam.Hidrodinamica_vs_normal_UP ,file="Hallmark_GENESET_HIDRODINAMICA_UP.Rdata")
write.csv(x=cam.Dieta_vs_normal_UP,"Hallmark_GENESET_HIDRODINAMICA_UP.csv")
save(cam.Hidrodinamica_vs_normal_Down,file="Hallmark_GENESET_HIDRODINAMICA_Down.Rdata")
write.csv(x=cam.Dieta_vs_normal_Down,"Hallmark_GENESET_HIDRODINAMICA_Down.csv")

#Transformados

idx <- ids2indices(Hs.H,id=Transformados_vs_normal_ortolog$human_entrez)
cam.Transformadas_vs_normal <- camera(v,idx,design,contrast=cont.matrix[,4])
#Ordenado
cam.Transformadas_vs_normal[order(cam.Transformadas_vs_normal$PValue),]
head(cam.Transformadas_vs_normal_UP,5)
#Up-regulados
cam.Transformadas_vs_normal_UP <- subset(x=cam.Transformadas_vs_normal,subset=Direction=="Up")
cam.Transformadas_vs_normal_UP[order(cam.Transformadas_vs_normal_UP$PValue),]
head(cam.Transformadas_vs_normal_UP,5)
#Down-regulados
cam.Transformadas_vs_normal_Down <- subset(x=cam.Transformadas_vs_normal,subset=Direction=="Down")
cam.Transformadas_vs_normal_Down[order(cam.Transformadas_vs_normal_Down$PValue),]
head(cam.Transformadas_vs_normal_Down,5)

            
save(cam.Transformadas_vs_normal,file="Hallmark_GENESET_TRANSFORMADAS.Rdata")
write.csv(x=cam.Transformadas_vs_normal,"Hallmark_GENESET_TRANSFORMADAS.csv")
save(cam.Transformadas_vs_normal_UP,file="Hallmark_GENESET_TRANSFORMADAS_UP.Rdata")
write.csv(x=cam.Transformadas_vs_normal_UP,"Hallmark_GENESET_TRANSFORMADAS_UP.csv")
save(cam.Transformadas_vs_normal_Down,file="Hallmark_GENESET_TRANSFORMADAS_Down.Rdata")
write.csv(x=cam.Transformadas_vs_normal_Down,"Hallmark_GENESET_TRANSFORMADAS_Down.csv")

#Dieta Vs hidrodinamica

idx <- ids2indices(Hs.H,id=DietaVsHidro_ortolog$human_entrez)
cam.DietaVsHidro <- camera(v,idx,design,contrast=cont.matrix[,5])
#Ordenarlos
cam.DietaVsHidro[order(cam.DietaVsHidro$PValue),]
head(cam.DietaVsHidro_UP,5)
#Up-regulados
cam.DietaVsHidro_UP <- subset(x=cam.DietaVsHidro,subset=Direction=="Up")
cam.DietaVsHidro_UP [order(cam.DietaVsHidro_UP $PValue),]
head(cam.DietaVsHidro_UP,5)
#Down-regulados
cam.DietaVsHidro_Down <- subset(x=cam.DietaVsHidro,subset=Direction=="Down")
cam.DietaVsHidro_Down[order(cam.DietaVsHidro_Down$PValue),]
head(cam.DietaVsHidro_Down,5)
            
save(cam.DietaVsHidro,file="Hallmark_GENESET_DietVsHidro.Rdata")
write.csv(x=cam.DietaVsHidro,"Hallmark_GENESET_DietVsHidro.csv")
save(cam.DietaVsHidro_UP,file="Hallmark_GENESET_DietVsHidro_UP.Rdata")
write.csv(x=cam.DietaVsHidro_UP,"Hallmark_GENESET_DietVsHidro_UP.csv")
save(cam.DietaVsHidro_Down,file="Hallmark_GENESET_DietVsHidro_Down.Rdata")
write.csv(x=cam.DietaVsHidro_Down,"Hallmark_GENESET_DietVsHidro_Down.csv")




