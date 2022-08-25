#Trabajo Fin de Máster - Máster en Métodos Computacionales en ciencias
#Autor: Oier Azurmendi Senar
#Tutores_ Silvestre Vicent Cambra y Mikel Hernaez
#Curso 2022-2023

#Objetivo: Construir el índice de los transcritos a partir del genoma de ratones GRCm39.

######################################################################################################


mkdir kallisto_getting_started

cd kallisto_getting_started

wget http://ftp.ensembl.org/pub/release-107/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz

gunzip Mus_musculus.GRCm39.cdna.all.fa.gz

kallisto index -i transcripts.idx Mus_musculus.GRCm39.cdna.all.fa.gz 
