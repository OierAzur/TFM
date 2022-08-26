library(babelgene)

#CCARAG vs NORMAL
CCARag_vs_normal_mvg_up <- CCARag_vs_normal[ order(CCARag_vs_normal$logFC,decreasing=T), ][1:200,1]
ortolog_CCARag_up <- orthologs(genes=CCARag_vs_normal_mvg_up,species="Mus musculus",human=F)
CCARagUP <- ortolog_CCARag_up$human_symbol
write.table(CCARagUP ,file="CCARagUP.txt",row.names = F,col.names = "F",quote = F)
CCARag_vs_normal_mvg_down <- CCARag_vs_normal[ order(CCARag_vs_normal$logFC), ][1:200,1]
ortolog_CCARag_down <- orthologs(genes=CCARag_vs_normal_mvg_down,species="Mus musculus",human=F)
CCARagDOWN <- ortolog_CCARag_down$human_symbol
write.table(CCARagDOWN ,file="CCARagDOWN.txt",row.names = F,col.names = "F",quote = F)

#Dieta VS Normal
Dieta_vs_normal_mvg_up <- Dieta_vs_normal[ order(Dieta_vs_normal$logFC,decreasing=T), ][1:200,1]
Dieta_vs_normal_mvg_down <- Dieta_vs_normal[ order(Dieta_vs_normal$logFC), ][1:200,1]
ortolog_Dieta_up <- orthologs(genes=Dieta_vs_normal_mvg_up,species="Mus musculus",human=F)
DietaUP <- ortolog_Dieta_up$human_symbol
write.table(DietaUP ,file="DietaUP.txt",row.names = F,col.names = "F",quote = F)
ortolog_Dieta_down <- orthologs(genes=Dieta_vs_normal_mvg_down,species="Mus musculus",human=F)
DietaDOWN <- ortolog_Dieta_down$human_symbol
write.table(DietaDOWN ,file="DietaDOWN.txt",row.names = F,col.names = "F",quote = F)

#Hidrodinamica vs normal

Hidrodinamica_vs_normal_mvg_up <- Hidrodinamica_vs_normal[ order(Hidrodinamica_vs_normal$logFC,decreasing=T), ][1:200,1]
Hidrodinamica_vs_normal_mvg_down <- Hidrodinamica_vs_normal[ order(Hidrodinamica_vs_normal$logFC), ][1:200,1]
ortolog_Hidrodinamica_up <- orthologs(genes=Hidrodinamica_vs_normal_mvg_up,species="Mus musculus",human=F)
HidrodinamicaUP <- ortolog_Hidrodinamica_up$human_symbol
write.table(HidrodinamicaUP ,file="HidrodinamicaUP.txt",row.names = F,col.names = "F",quote = F)
ortolog_Hidrodinamica_down <- orthologs(genes=Hidrodinamica_vs_normal_mvg_down,species="Mus musculus",human=F)
HidrodinamicaDOWN <- ortolog_Hidrodinamica_down$human_symbol
write.table(HidrodinamicaDOWN ,file="HidrodinamicaDOWN.txt",row.names = F,col.names = "F",quote = F)

#Celulas Transformadas vs normal
Transformadas_vs_normal
Transformadas_vs_normal_mvg_up <- Transformadas_vs_normal[ order(Transformadas_vs_normal$logFC,decreasing=T), ][1:200,1]
Transformadas_vs_normal_mvg_down <-Transformadas_vs_normal[ order(Transformadas_vs_normal$logFC), ][1:200,1]
ortolog_Transformadas_up <- orthologs(genes=Transformadas_vs_normal_mvg_up,species="Mus musculus",human=F)
TransformadasUP <- ortolog_Transformadas_up$human_symbol
write.table(TransformadasUP ,file="TransformadasUP.txt",row.names = F,col.names = "F",quote = F)
ortolog_Transformadas_down <- orthologs(genes=Transformadas_vs_normal_mvg_down,species="Mus musculus",human=F)
TransformadasDOWN <- ortolog_Transformadas_down$human_symbol
write.table(TransformadasDOWN ,file="TransformadasDOWN.txt",row.names = F,col.names = "F",quote = F)

#Dieta vs Hidrodinamica
Dieta_vs_Hidrodinamica
Dieta_vs_Hidrodinamica_mvg_up <- Dieta_vs_Hidrodinamica[ order(Dieta_vs_Hidrodinamica$logFC,decreasing=T), ][1:200,1]
Dieta_vs_Hidrodinamica_mvg_down <-Dieta_vs_Hidrodinamica[ order(Dieta_vs_Hidrodinamica$logFC), ][1:200,1]
ortolog_Dieta_vs_Hidrodinamica_up <- orthologs(genes=Dieta_vs_Hidrodinamica_mvg_up,species="Mus musculus",human=F)
DIETAvsHidrodinamicaUP <- ortolog_Dieta_vs_Hidrodinamica_up$human_symbol
write.table(DIETAvsHidrodinamicaUP  ,file="DIETAvsHidrodinamicaUP .txt",row.names = F,col.names = "F",quote = F)
ortolog_Dieta_vs_Hidrodinamica_down<- orthologs(genes=Dieta_vs_Hidrodinamica_mvg_down,species="Mus musculus",human=F)
DIETAvsHidrodinamicaDOWN <- ortolog_Dieta_vs_Hidrodinamica_down$human_symbol
write.table(DIETAvsHidrodinamicaDOWN  ,file="DIETAvsHidrodinamicaDOWN .txt",row.names = F,col.names = "F",quote = F)
