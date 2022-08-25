#Trabajo Fin de Máster - Máster en Métodos Computacionales en ciencias
#Autor: Oier Azurmendi Senar
#Tutores_ Silvestre Vicent Cambra y Mikel Hernaez
#Curso 2021-2022

#Objetivo: Primera observación preliminar de las muestras y condiciones, y análisis de expresión 
#diferencial de los isoformas de los genes. 

######################################################################################################
setwd("Documents/TFM/Kallisto_complete_Results/")
library("gridExtra")
library("cowplot")
library(rhdf5)
library("sleuth")
library("biomaRt")
library(tidyverse)

#Cargamos los datos
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

#Incluir los nombres de los genes en Análisis a nivel de transcrito
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "mmusculus_gene_ensembl",
                         host = 'ensembl.org')

t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                     "external_gene_name"), mart = mart)

t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
save(t2g,file="Biomart_DGI.Rdata")

#Sleuth a nivel de transcrito:
so <- sleuth_prep(s2c, 
                  full_model = ~condition, 
                  target_mapping = t2g, 
                  read_bootstrap_tpm = TRUE,
                  extra_bootstrap_summary = TRUE,
                  transformation_function = function(x) log2(x + 0.5))

#Se perfila el modelo acorde a las condiciones
so <- sleuth_fit(so)
so <- sleuth_fit(so, ~condition, 'full')
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
models(so)

#PCA
plot_pca(so, 
         color_by = 'condition',
         text_labels = TRUE)

#Varianza
plot_pc_variance(so)

#Loadings
plot_loadings(so)

#Heatmap de las muestras
plot_sample_heatmap(so)
                                                     
#Analisis de los resultados
sleuth_results_so <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)

#Significativos
sig_transcripts <- sleuth_results_so %>% 
  filter(qval < 0.05)
head(sig_transcripts, 20)
write.csv(sig_transcripts="DGI_lrt_passed.csv")

Por ejemplo miro el primero según las condiciones
plot_bootstrap(so, "ENSMUST00000053361.12", units = "est_counts", color_by = "condition")
plot_bootstrap(so, "ENSMUST00000078357.5", units = "est_counts", color_by = "condition")
plot_bootstrap(so, "ENSMUST00000193675.2", units = "est_counts", color_by = "condition")

#Transcritos más significados
plot_transcript_heatmap(so, 
                        transcripts = sig_transcripts$target_id[1:20])



#Se observan las densidades de cada grupo
plot_group_density(so, use_filtered = TRUE, units = "est_counts",
                   trans = "log", grouping = setdiff(colnames(so$sample_to_covariates),
                                                     "sample"), offset = 1)

plot_scatter(so)


#Obtención de la matriz
matriz_sleuth <- sleuth_to_matrix(so_g,'obs_norm','tpm')
df_sleuth <- data.frame(matriz_sleuth)
df_sleuth <- tibble::rownames_to_column(df_sleuth)

#Obtener HTML
sleuth_live(so_g)
