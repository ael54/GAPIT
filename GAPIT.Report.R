`GAPIT.Report` <-
function(name.of.trait=NULL,GWAS=NULL,pred=pred,tvalue=NULL,stderr=NULL,Vp=1,
DPP=100000,cutOff=.01,threshold.output=.01,MAF=NULL,seqQTN=NULL){
#Object: Out put plots and tables
#Input: GWAS,name.of.trait, DPP 
#Requirement: None
#Output: Graphs and tables
#Authors: Zhiwu Zhang
# Date  start: April 2, 2013
# Last update: April 2, 2013
##############################################################################################
#print("GAPIT.Report Started")
#print(seqQTN)
#Manhattan Plots
#print("Manhattan plot (Genomewise)..." )
GAPIT.Manhattan(GI.MP = GWAS[,2:4], name.of.trait = name.of.trait, DPP=DPP, plot.type = "Genomewise",cutOff=cutOff,seqQTN=seqQTN)
#print("Manhattan plot (Chromosomewise)..." )
GAPIT.Manhattan(GI.MP = GWAS[,2:4], name.of.trait = name.of.trait, DPP=DPP, plot.type = "Chromosomewise",cutOff=cutOff)

#QQ plots
#print("QQ plotting..." )
GAPIT.QQ(P.values = GWAS[,4], name.of.trait = name.of.trait,DPP=DPP)

#Association Table
#print("Create association table..." )
index=1:nrow(GWAS)
if(threshold.output<1)index=which(GWAS[,4]<threshold.output)
write.table(GWAS[index,], paste("GAPIT.", name.of.trait, ".GWAS.Results.csv", sep = ""), quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)

#Prediction
#print("Create prediction table..." )
if(!is.null(pred)) write.table(pred, paste("GAPIT.", name.of.trait, ".Pred.csv", sep = ""), quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)

#ROC
#print("Creating ROC table and plot" )
myROC=GAPIT.ROC(t=tvalue,se=stderr,Vp=Vp,trait=name.of.trait)

#MAF
#print("Creating MAF table and plot" )
myMAF=GAPIT.MAF(MAF=MAF,P=GWAS[,4],E=NULL,trait=name.of.trait,threshold.output=threshold.output)

#print("Report accomplished" )
}#The function GAPIT.Report ends here
