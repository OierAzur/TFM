#Trabajo Fin de Máster - Máster en Métodos Computacionales en ciencias
#Autor: Oier Azurmendi Senar
#Tutores_ Silvestre Vicent Cambra y Mikel Hernaez
#Curso 2021-2022

#Objetivo: Análisis de supervivencia a partir de datos clínicos del TCGA-CHOL.

######################################################################################################
setwd("Documents/TFM/Surv/")


library("org.Hs.eg.db")
library("TCGAbiolinks")
library("TCGAbiolinksGUI.data")
library("limma")
library("edgeR")
library("glmnet")
library("factoextra")
library("SummarizedExperiment")
library("gplots")
library("survival")
library("survminer")
library("RColorBrewer")
library("gProfileR")
library("genefilter")


#----------TCGA data-----------
GDCprojects = getGDCprojects()

head(GDCprojects[c("project_id", "name")])

#Resumen de los datos de Coangiocarcinoma
TCGAbiolinks:::getProjectSummary("TCGA-CHOL")

# Datos de expresión
query_TCGA = GDCquery(
  project = "TCGA-CHOL",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  data.type = "Gene Expression Quantification",
  experimental.strategy = "RNA-Seq",
  sample.type = c("Primary Tumor", "Solid Tissue Normal"),
  access = "open")

# Metadata de los datos
metadata <- query_TCGA[[1]][[1]]

#En windows se pueden descargar los datos directamente.
GDCdownload(query = query_TCGA)
tcga_data = GDCprepare(query_TCGA)

#Desde Mac el comando GDCprepare no funciona. Se han descargado en windows, descargado y cargado aquí:
load("/Users/oierazur/Downloads/datosrurv.Rdata")

#Parámetros de los datos clínicos.
colnames(colData(tcga_data))
table(tcga_data@colData$vital_status)
table(tcga_data@colData$definition)
table(tcga_data@colData$tissue_or_organ_of_origin)
table(tcga_data@colData$gender)
table(tcga_data@colData$race)

#Matrices de expresión Genica
dim(assay(tcga_data))  
head(rowData(tcga_data)) 


limma_pipeline = function(
  datos_tcga,
  variable_interes,
  grupo_referencia=NULL){
  
  factor_diseno = colData(datos_tcga)[, variable_interes, drop=T]
  
  grupo = factor(factor_diseno)
  if(!is.null(grupo_referencia)){grupo = relevel(grupo, ref=grupo_referencia)}
  
  diseno = model.matrix(~ grupo)
  
  dge = DGEList(counts=assay(datos_tcga),
                samples=colData(datos_tcga),
                genes=as.data.frame(rowData(datos_tcga)))
  
  # filtrado
  keep = filterByExpr(dge,diseno)
  dge = dge[keep,,keep.lib.sizes=FALSE]
  rm(keep)
  
  # Normalización 
  dge = calcNormFactors(dge)
  v = voom(dge, diseno, plot=TRUE)
  
  # Hacer los contrastes
  fit = lmFit(v, diseno)
  fit = eBayes(fit)
  
  # Obtencion de tablas
  topGenes = topTable(fit, coef=ncol(diseno),number=20000, sort.by="p")
  
  return(
    list(
      voomObj=v, # datos normalizados
      fit=fit, # modelo lineal y estadistivas
      topGenes=topGenes # 100 genes mas DGE 
    )
  )
}

limma_res = limma_pipeline(
  datos_tcga=tcga_data,
  variable_interes="definition",
  grupo_referencia="Solid Tissue Normal"
)

#Datos Clínicos

clinical = tcga_data@colData

dim(clinical)

#Obtenemos las variables de interes
clin_df = clinical[clinical$definition == c("Primary solid Tumor","Solid Tissue Normal"),
                   c("patient",
                     "vital_status",
                     "days_to_death",
                     "days_to_last_follow_up",
                     "gender",
                     "ajcc_pathologic_stage",
                     "age_at_diagnosis")]

# Creamos un valor booleano: True para pacientes fallecidos
#  y False para pacientes vivos
clin_df$deceased = clin_df$vital_status == "Dead"

clin_df$overall_survival = ifelse(clin_df$deceased,
                                  clin_df$days_to_death,
                                  clin_df$days_to_last_follow_up)

head(clin_df)

#Expresión génica y supervivencia
expr_df = limma_res$topGenes

print(expr_df[,1 ])

# Obtención de las matrices 
d_mat = as.matrix(t(limma_res$voomObj$E))
print(dim(d_mat))
# En factor
d_resp = as.factor(limma_res$voomObj$targets$definition)

#Introducimos el nombre del gen en interes en cada estudio
gene_interes <- "KLF6"
index_gene <- which(expr_df$gene_name==gene_interes)
gene_id = expr_df[index_gene, "gene_id"]
gene_name = expr_df[index_gene, "gene_name"]

expr_diseased = d_mat[rownames(clin_df), gene_id]
expr_healthy = d_mat[setdiff(rownames(d_mat), rownames(clin_df)), gene_id]

#Comparamos los niveles de expresión en sanos vs enfermedad
boxplot(expr_diseased, expr_healthy,
        names=c("Diseased", "Healthy"), main="Distribution of gene expression")

# Valores de expresión del gen seleccionado
clin_df$gene_value = d_mat[rownames(clin_df), gene_id]

#La mediana
median_value = median(clin_df$gene_value)
print(median_value)

# División de los pacientes acorde a la expresión del gen
clin_df$gene = ifelse(clin_df$gene_value >= median_value, "UP", "DOWN")

# Modelo de supervivencia
fit = survfit(Surv(overall_survival, deceased)~gene, data=clin_df)

# P-valor del modelo de supervivencia
pval = surv_pvalue(fit, data=clin_df)$pval
print(pval)

#Obtenemos las tablas
surv_summary(fit,data=clin_df)

# Kaplan-Meier plot
ggsurvplot(fit, data=clin_df, legend.title="Grupo",
           pval=T, risk.table=T, title=paste(gene_name),legend="none")

#Log-rank test
surv_diff <- survdiff(Surv(overall_survival, deceased) ~gene, data=clin_df)
surv_diff

#Multivariate con el gen de interes
ggsurvplot(survfit(Sur~clin_df$gene+estado), data=clin_df, legend.title="Grupo",
           pval=T, risk.table=T, title=paste(gene_name),xlab="Time in days",
           break.time.by = 365,legend="none",
           legend.labs=c("Stage I down regulated","Stage II down regulated",
                         "Stage III down regulated",
                         "Stage IVB down regulated","Stage I up regulated",
                         "Stage II up regulated",
                         "Stage IVB up regulated"),ggtheme = theme_gray())


##-----------Con valores de todos los genes-----------##

####Obtener tablas de valores para poder ordenar
#Obtener el d_mat solo con los pacientes que tenemos en clin_df 
d_mat = as.matrix(limma_res$voomObj$E)

#Recortar nombre de paciente en d_mat
rownames(d_mat) <- lapply(rownames(d_mat),sub,pattern="\\.\\d+$",replacement="")

## Obtener Cox regression para cada gen
Res = data.frame(Gene = rownames(d_mat), LHR=NA, PV = 1, FDR=1)
rownames(Res) = rownames(d_mat)
ig=1

keep <-which(clin_df@rownames%in%colnames(d_mat))
d_mat<- d_mat[,keep]
Sur = Surv(time = clin_df$overall_survival, event =clin_df$deceased)
for (ig in 1:nrow(d_mat)){
  mod = coxph(Sur ~ d_mat[ig,])
  Res$LHR[ig] = summary(mod)$coefficients[,"coef"]
  Res$PV[ig] = summary(mod)$sctest["pvalue"]
  
  ## write message every 1000 genes
  if(ig%%1000 == 0){
    print(ig)
    flush.console()
  }
}

Res$FDR = p.adjust(Res$PV,method="fdr")

str(Res)

#Obtenemos la lista de genes con un p-valor menor a 0.05
keep<- Res$PV<0.05
Res <- Res[keep,]

OrgDB <- org.Hs.eg.db
Res$Gene <-  mapIds(OrgDB,keys=Res$Gene,
                          column="SYMBOL",
                          keytype="ENSEMBL",
                          multiVals="first")
#Los ordenamos
genes_analizar <- genes_analizar[ order(genes_analizar$PV ,decreasing=T), ]
isig = Res$FDR<0.1
Res[isig,-1]

# Combinamos las expresiones de todos los genes
score = double(ncol(d_mat))

ip=1
for (ip in 1:ncol(d_mat)){
  score[ip] = sum(Res$LHR[isig] *d_mat[isig,ip])
}

#Vemos las relaciones entre variantes
mod <- coxph(Sur~score)
summary(mod)

#Kaplan - Meier plot
#Al ser análisis multivariados, con cada factor:
#Sexo:
expresion= c("low","high")[1+as.integer(score>median(score))]
expresion = factor(expresion,levels=c("low","high"))
sexo=as.factor(clin_df$gender)
ggsurvplot(survfit(Sur ~ expresion+sexo), data=clin_df,legend="none",
           legend.title="Grupo",
           pval=T, risk.table=T, title="Sexo")
surv_summary(survfit(Sur ~ expresion+sexo),data=clin_df)

#Con el estado del tumor:
estado=as.factor(clin_df$ajcc_pathologic_stage)
ggsurvplot(survfit(Sur ~ expresion+estado), data=clin_df,legend="none",
           legend.title="Grupo",
           pval=T, risk.table=T, title="Stage")
surv_summary(survfit(Sur ~ expresion+estado),data=clin_df)

#Con la edad:
edad = c("<65 años",">65 años")[1+as.integer(clin_df$age_at_diagnosis/365<65)]
edad = factor(edad,levels=c("<65 años",">65 años"))
ggsurvplot(survfit(Sur ~ expresion+edad), data=clin_df, legend.title="Grupo",legend="none",
           pval=T, risk.table=T, title="Edad")
surv_summary(survfit(Sur ~ expresion+edad),data=clin_df)




           
