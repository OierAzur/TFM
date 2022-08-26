
# Load packages
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
#https://docs.gdc.cancer.gov/Data_Portal/Users_Guide/Projects/

TCGAbiolinks:::getProjectSummary("TCGA-CHOL")
##the summary shows that TCGA provide data for 500 patients including clinical, expression, DNA methylation, and genotyping data

##for simplicity, we ignore the small class of metastasis, therefore, we redo the query
query_TCGA = GDCquery(
  project = "TCGA-CHOL",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  data.type = "Gene Expression Quantification",
  experimental.strategy = "RNA-Seq",
  sample.type = c("Primary Tumor", "Solid Tissue Normal"),
  access = "open")

# Get metadata matrix
metadata <- query_TCGA[[1]][[1]]
##Next, we need to download the files from the query
GDCdownload(query = query_TCGA)
tcga_data = GDCprepare(query_TCGA)

#Desde Mac
load("/Users/oierazur/Downloads/datosrurv.Rdata")

colnames(colData(tcga_data))
table(tcga_data@colData$vital_status)
table(tcga_data@colData$definition)
table(tcga_data@colData$tissue_or_organ_of_origin)
table(tcga_data@colData$gender)
table(tcga_data@colData$race)

#Gene expression matrices
dim(assay(tcga_data))  

head(rowData(tcga_data)) 

limma_pipeline = function(
  tcga_data,
  condition_variable,
  reference_group=NULL){
  
  design_factor = colData(tcga_data)[, condition_variable, drop=T]
  
  group = factor(design_factor)
  if(!is.null(reference_group)){group = relevel(group, ref=reference_group)}
  
  design = model.matrix(~ group)
  
  dge = DGEList(counts=assay(tcga_data),
                samples=colData(tcga_data),
                genes=as.data.frame(rowData(tcga_data)))
  
  # filtering
  keep = filterByExpr(dge,design)
  dge = dge[keep,,keep.lib.sizes=FALSE]
  rm(keep)
  
  # Normalization (TMM followed by voom)
  dge = calcNormFactors(dge)
  v = voom(dge, design, plot=TRUE)
  
  # Fit model to data given design
  fit = lmFit(v, design)
  fit = eBayes(fit)
  
  # Show top genes
  topGenes = topTable(fit, coef=ncol(design),number=20000, sort.by="p")
  
  return(
    list(
      voomObj=v, # normalized data
      fit=fit, # fitted model and statistics
      topGenes=topGenes # the 100 most differentially expressed genes
    )
  )
}

limma_res = limma_pipeline(
  tcga_data=tcga_data,
  condition_variable="definition",
  reference_group="Solid Tissue Normal"
)
###Clinical Data

# extract clinical data
clinical = tcga_data@colData

dim(clinical)
clin_df = clinical[clinical$definition == c("Primary solid Tumor","Solid Tissue Normal"),
                   c("patient",
                     "vital_status",
                     "days_to_death",
                     "days_to_last_follow_up",
                     "gender",
                     "ajcc_pathologic_stage",
                     "age_at_diagnosis")]

# create a new boolean variable that has TRUE for dead patients
# and FALSE for live patients
clin_df$deceased = clin_df$vital_status == "Dead"

clin_df$overall_survival = ifelse(clin_df$deceased,
                                  clin_df$days_to_death,
                                  clin_df$days_to_last_follow_up)

# show first 10 samples
head(clin_df)


#Gene expression and Survival
expr_df = limma_res$topGenes

print(expr_df[,1 ])

# visualize the gene expression distribution on the diseased samples (in black)
# versus the healthy samples (in red)

d_mat = as.matrix(t(limma_res$voomObj$E))
print(dim(d_mat))
# As before, we want this to be a factor
d_resp = as.factor(limma_res$voomObj$targets$definition)

#Introducimos el nombre del gen en interes
gene_interes <- "KLF6"
index_gene <- which(expr_df$gene_name==gene_interes)
# also get the common gene name of the first row
gene_id = expr_df[index_gene, "gene_id"]
gene_name = expr_df[index_gene, "gene_name"]

# visualize the gene expression distribution on the diseased samples (in black)
# versus the healthy samples (in red)
expr_diseased = d_mat[rownames(clin_df), gene_id]
expr_healthy = d_mat[setdiff(rownames(d_mat), rownames(clin_df)), gene_id]

boxplot(expr_diseased, expr_healthy,
        names=c("Diseased", "Healthy"), main="Distribution of gene expression")

# get the expression values for the selected gene
clin_df$gene_value = d_mat[rownames(clin_df), gene_id]

median_value = median(clin_df$gene_value)
print(median_value)

# divide patients in two groups, up and down regulated.
# if the patient expression is greater or equal to them median we put it
# among the "up-regulated", otherwise among the "down-regulated"
clin_df$gene = ifelse(clin_df$gene_value >= median_value, "UP", "DOWN")

# we can fit a survival model, like we did in the previous section
fit = survfit(Surv(overall_survival, deceased)~gene, data=clin_df)

# we can extract the survival p-value and print it
pval = surv_pvalue(fit, data=clin_df)$pval
print(pval)

#Obtenemos las tablas
surv_summary(fit,data=clin_df)

# and finally, we produce a Kaplan-Meier plot
ggsurvplot(fit, data=clin_df, legend.title="Grupo",
           pval=T, risk.table=T, title=paste(gene_name),legend="none")

#Log-rank test
surv_diff <- survdiff(Surv(overall_survival, deceased) ~gene, data=clin_df)
surv_diff

#Multivariate con un gen en concreto
ggsurvplot(survfit(Sur~clin_df$gene+estado), data=clin_df, legend.title="Grupo",
           pval=T, risk.table=T, title=paste(gene_name),xlab="Time in days",
           break.time.by = 365,legend="none",ggtheme = theme_gray())
ggsurvplot(survfit(Sur~clin_df$gene+estado), data=clin_df, legend.title="Grupo",
           pval=T, risk.table=T, title=paste(gene_name),xlab="Time in days",
           break.time.by = 365,legend="none",
           legend.labs=c("Stage I down regulated","Stage II down regulated",
                         "Stage III down regulated",
                         "Stage IVB down regulated","Stage I up regulated",
                         "Stage II up regulated",
                         "Stage IVB up regulated"),ggtheme = theme_gray())


####Obtener tablas de valores para poder ordenar

#Obtener el d_mat solo con los pacientes que tenemos en clin_df 
#Recortar nombre de paciente en d_mat
d_mat = as.matrix(limma_res$voomObj$E)
rownames(d_mat) <- lapply(rownames(d_mat),sub,pattern="\\.\\d+$",replacement="")

## build container and Cox regression for each gene
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
#keep<- Res$PV<0.05
#Res <- Res[keep,]
library(org.Hs.eg.db)
OrgDB <- org.Hs.eg.db
Res$Gene <-  mapIds(OrgDB,keys=Res$Gene,
                          column="SYMBOL",
                          keytype="ENSEMBL",
                          multiVals="first")
genes_analizar <- genes_analizar[ order(genes_analizar$PV ,decreasing=T), ]
isig = Res$FDR<0.1
Res[isig,-1]


## lets combine all the genes expresion######
score = double(ncol(d_mat))

ip=1
for (ip in 1:ncol(d_mat)){
  score[ip] = sum(Res$LHR[isig] *d_mat[isig,ip])
}

mod <- coxph(Sur~score)
summary(mod)


## make a KM plot
#score -> factor
expresion= c("low","high")[1+as.integer(score>median(score))]
expresion = factor(expresion,levels=c("low","high"))
sexo=as.factor(clin_df$gender)
ggsurvplot(survfit(Sur ~ expresion+sexo), data=clin_df,legend="none",
           legend.title="Grupo",
           pval=T, risk.table=T, title="Sexo")
surv_summary(survfit(Sur ~ expresion+sexo),data=clin_df)
#Con el estado del tumor
estado=as.factor(clin_df$ajcc_pathologic_stage)
ggsurvplot(survfit(Sur ~ expresion+estado), data=clin_df,legend="none",
           legend.title="Grupo",
           pval=T, risk.table=T, title="Stage")
surv_summary(survfit(Sur ~ expresion+estado),data=clin_df)
#Con la edad
edad = c("<65 años",">65 años")[1+as.integer(clin_df$age_at_diagnosis/365<65)]
edad = factor(edad,levels=c("<65 años",">65 años"))
ggsurvplot(survfit(Sur ~ expresion+edad), data=clin_df, legend.title="Grupo",legend="none",
           pval=T, risk.table=T, title="Edad")
surv_summary(survfit(Sur ~ expresion+edad),data=clin_df)




           
