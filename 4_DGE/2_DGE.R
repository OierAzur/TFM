#Trabajo Fin de Máster - Máster en Métodos Computacionales en ciencias
#Autor: Oier Azurmendi Senar
#Tutores_ Silvestre Vicent Cambra y Mikel Hernaez
#Curso 2021-2022

#Objetivo:Análisis de la expresión diferencial de los genes

######################################################################################################
setwd("Documents/TFM/DGE/")
library(limma)
library(edgeR) 
library(AnnotationDbi)

#Le aplicamos una matriz de diseño y de contraste
group <- factor(s2c$condition)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
design

cont.matrix <- makeContrasts(CCARag2=CCA.Rag2_F1  - Colangiocitos.normales,Dieta=Dieta - Colangiocitos.normales,
                             Hidrodinamica=Hidrodinamica-Colangiocitos.normales,
                             Transformados=Colangiocitos.transformados - Colangiocitos.normales,
                             DietaVsHidro=Dieta-Hidrodinamica,levels=design)
cont.matrix

#Aplicamos voom
par(mfrow=c(1,1))
v <- voom(y,design,plot = TRUE)
names(v)

#Test para expresión diferencial mediante modelos lineares
fit <- lmFit(v,design)
names(fit)
head(coef(fit))

#Creamos matriz de conteo
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
dim(fit.cont)
summa.fit <- decideTests(fit.cont,p.value  = 0.05,lfc=1)
summary(summa.fit)

#Guardar tabla
save(fit.cont,file="fit_DGE.Rdata")
save(v,file="voom_result.Rdata")

#Gráficas
#Venn 
vennDiagram(summa.fit,include=c("up", "down"),
counts.col=c("red", "blue"),
circle.col = c("red", "blue", "green3","yellow"))

#Diferencía de medianas de los datos de expresión
plotMD(fit.cont, column=1, status=summa.fit[,1], main=colnames(fit.cont)[1])
plotMD(fit.cont, column=2, status=summa.fit[,2], main=colnames(fit.cont)[2])
plotMD(fit.cont, column=3, status=summa.fit[,3], main=colnames(fit.cont)[3])
plotMD(fit.cont, column=4, status=summa.fit[,4], main=colnames(fit.cont)[4])
plotMD(fit.cont, column=4, status=summa.fit[,5], main=colnames(fit.cont)[5])
