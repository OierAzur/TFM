#Trabajo Fin de Máster - Máster en Métodos Computacionales en ciencias
#Autor: Oier Azurmendi Senar
#Tutores_ Silvestre Vicent Cambra y Mikel Hernaez
#Curso 2022-2023

#Objetivo: Ejecutar el algorítmo de cuantificación de los transcritos mediante el pseudo alineamiento.

######################################################################################################

mkdir kallisto_complete_results

cd kallisto_complete results

#Muestras de Colangiocitos Normales
kallisto quant -i transcripts.idx -o KPE1_S13 -b 100 --single -l 51 -s 20 KPE1_S13_R1_001.fastq 
kallisto quant -i transcripts.idx -o KPE2_S14 -b 100 --single -l 51 -s 20 KPE2_S14_R1_001.fastq 
kallisto quant -i transcripts.idx -o KPE3_S15 -b 100 --single -l 51 -s 20 KPE3_S15_R1_001.fastq

#Muestras de Colangiocitos Transformados
kallisto quant -i transcripts.idx -o KPC1_S10 -b 100 --single -l 51 -s 20 KPC1_S10_R1_001.fastq 
kallisto quant -i transcripts.idx -o KPC2_S11 -b 100 --single -l 51 -s 20 KPC2_S11_R1_001.fastq 
kallisto quant -i transcripts.idx -o KPC3_S12 -b 100 --single -l 51 -s 20 KPC3_S12_R1_001.fastq 

#Muestras de CCA (Rag-F1)
kallisto quant -i transcripts.idx -o C_17741_S7 -b 100 --single -l 51 -s 20 C_17741_S7_R1_001.fastq 
kallisto quant -i transcripts.idx -o C_17742_S8 -b 100 --single -l 51 -s 20 C_17742_S8_R1_001.fastq 
kallisto quant -i transcripts.idx -o C_17743_S9 -b 100 --single -l 51 -s 20 C_17743_S9_R1_001.fastq

#Muestras de Hidrodinámica
kallisto quant -i transcripts.idx -o H_20791_S1 -b 100 --single -l 51 -s 20 H_20791_S1_R1_001.fastq 
kallisto quant -i transcripts.idx -o H_20792_S2 -b 100 --single -l 51 -s 20 H_20792_S2_R1_001.fastq 
kallisto quant -i transcripts.idx -o H_20793_S3 -b 100 --single -l 51 -s 20 H_20793_S3_R1_001.fastq

#Muestras de Dieta
kallisto quant -i transcripts.idx -o D_8161_S4 -b 100 --single -l 51 -s 20 D_8161_S4_R1_001.fastq 
kallisto quant -i transcripts.idx -o D_8161_S5 -b 100 --single -l 51 -s 20 D_8162_S5_R1_001.fastq 
kallisto quant -i transcripts.idx -o D_8161_S6 -b 100 --single -l 51 -s 20 D_8163_S6_R1_001.fastq





