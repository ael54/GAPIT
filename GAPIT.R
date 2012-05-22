`GAPIT` <-
function(Y=NULL,G=NULL,GD=NULL,GM=NULL,KI=NULL,Z=NULL,CV=NULL,CV.Inheritance=NULL,GP=NULL,GK=NULL,
                group.from=30 ,group.to=1000000,group.by=10,DPP=100000, 
                kinship.cluster="average", kinship.group='Mean',kinship.algorithm="VanRaden",                                                    
                bin.from=10000,bin.to=10000,bin.by=10000,inclosure.from=10,inclosure.to=10,inclosure.by=10,
                SNP.P3D=TRUE,SNP.effect="Add",SNP.impute="Middle",PCA.total=0, SNP.fraction = 1, seed = 123, BINS = 20,SNP.test=TRUE, 
                SNP.MAF=0,FDR.Rate = 1, SNP.FDR=1,SNP.permutation=FALSE,SNP.CV=NULL,SNP.robust="GLM",                             
                file.from=1, file.to=1, file.total=NULL, file.fragment = 99999,file.path=NULL, 
                file.G=NULL, file.Ext.G=NULL,file.GD=NULL, file.GM=NULL, file.Ext.GD=NULL,file.Ext.GM=NULL, 
                ngrid = 100, llim = -10, ulim = 10, esp = 1e-10,
                LD.chromosome=NULL,LD.location=NULL,LD.range=NULL,
                sangwich.top=NULL,sangwich.bottom=NULL,QC=TRUE,GTindex=NULL,LD=0.01,
                file.output=TRUE,cutOff=0.01, Model.selection = FALSE,output.numerical = FALSE,
                output.hapmap = FALSE, Create.indicator = FALSE,
				QTN=NULL, QTN.round=1,QTN.limit=0, QTN.update=TRUE, QTN.method="Penalty"){
#Object: To perform GWAS and GPS (Genomic Prediction/Selection)
#Designed by Zhiwu Zhang
#Writen by Alex Lipka, Feng Tian and Zhiwu Zhang
#Last update: October 24, 2011 
##############################################################################################
print("--------------------- Welcome to GAPIT ----------------------------")
echo=TRUE
GAPIT.Version=GAPIT.0000()

Timmer=GAPIT.Timmer(Infor="GAPIT")
Memory=GAPIT.Memory(Infor="GAPIT")

#Genotype processing and calculation Kin and PC
#First call to genotype to setup genotype data

myGenotype<-GAPIT.Genotype(G=G,GD=GD,GM=GM,KI=KI,kinship.algorithm=kinship.algorithm,PCA.total=PCA.total,SNP.fraction=SNP.fraction,SNP.test=SNP.test,
                file.path=file.path,file.from=file.from, file.to=file.to, file.total=file.total, file.fragment = file.fragment, file.G=file.G, 
                file.Ext.G=file.Ext.G,file.GD=file.GD, file.GM=file.GM, file.Ext.GD=file.Ext.GD,file.Ext.GM=file.Ext.GM,
                SNP.MAF=SNP.MAF,FDR.Rate = FDR.Rate,SNP.FDR=SNP.FDR,SNP.effect=SNP.effect,SNP.impute=SNP.impute,
                LD.chromosome=LD.chromosome,LD.location=LD.location,LD.range=LD.range,
                GP=GP,GK=GK,bin.size=NULL,inclosure.size=NULL, Timmer = Timmer,Memory=Memory,
                sangwich.top=sangwich.top,sangwich.bottom=sangwich.bottom,GTindex=NULL,file.output=file.output, Create.indicator = Create.indicator)

Timmer=myGenotype$Timmer
Memory=myGenotype$Memory

Timmer=GAPIT.Timmer(Timmer=Timmer,Infor="Genotype for all")
Memory=GAPIT.Memory(Memory=Memory,Infor="Genotype for all")


KI=myGenotype$KI
PC=myGenotype$PC

genoFormat=myGenotype$genoFormat
hasGenotype=myGenotype$hasGenotype
byFile=myGenotype$byFile
fullGD=myGenotype$fullGD
GD=myGenotype$GD
GI=myGenotype$GI
GT=myGenotype$GT
G=myGenotype$G

rownames(GD)=GT
colnames(GD)=GI[,1]

if(output.numerical) write.table(GD,  "GAPIT.Genotype.Numerical.txt", quote = FALSE, sep = "\t", row.names = TRUE,col.names = NA)
if(output.hapmap) write.table(myGenotype$G,  "GAPIT.Genotype.hmp.txt", quote = FALSE, sep = "\t", row.names = FALSE,col.names = FALSE)

#In case of null Y and null GP, return genotype only  
if(is.null(Y) & is.null(GP)) return (list(GWAS=NULL,GPS=NULL,Pred=NULL,compression=NULL,kinship.optimum=NULL,kinship=myGenotype$KI,PCA=myGenotype$PC,GD=data.frame(cbind(as.data.frame(GT),as.data.frame(GD))),GI=GI,G=myGenotype$G))

#In case of null Y, return genotype only          
if(is.null(Y)) return (list(GWAS=NULL,GPS=NULL,Pred=NULL,compression=NULL,kinship.optimum=NULL,kinship=myGenotype$KI,PCA=myGenotype$PC,GD=data.frame(cbind(as.date.frame(GT),as.data.frame(GD))),Gi=GI,G=myGenotype$G))

rm(myGenotype)
gc()

print("--------------------Processing traits----------------------------------")
if(!is.null(Y)){
print("Phenotype provided!")
if(ncol(Y)<2)  stop ("Phenotype should have taxa name and one trait at least. Please correct phenotype file!")

for (trait in 2: ncol(Y))  {
print(paste("Processing trait: ",colnames(Y)[trait],sep=""))
gapitMain <- GAPIT.Main(Y=Y[,c(1,trait)],G=G,GD=GD,GM=GM,KI=KI,Z=Z,CV=CV,CV.Inheritance=CV.Inheritance,GP=GP,GK=GK,SNP.P3D=SNP.P3D,kinship.algorithm=kinship.algorithm,
                      bin.from=bin.from,bin.to=bin.to,bin.by=bin.by,inclosure.from=inclosure.from,inclosure.to=inclosure.to,inclosure.by=inclosure.by,
				              group.from=group.from,group.to=group.to,group.by=group.by,kinship.cluster=kinship.cluster,kinship.group=kinship.group,name.of.trait = colnames(Y)[trait],
                        file.path=file.path,file.from=file.from, file.to=file.to, file.total=file.total, file.fragment = file.fragment, file.G=file.G,file.Ext.G=file.Ext.G,file.GD=file.GD, file.GM=file.GM, file.Ext.GD=file.Ext.GD,file.Ext.GM=file.Ext.GM, 
                        SNP.MAF= SNP.MAF,FDR.Rate = FDR.Rate,SNP.FDR=SNP.FDR,SNP.effect=SNP.effect,SNP.impute=SNP.impute,PCA.total=PCA.total,GAPIT.Version=GAPIT.Version,
                        GT=GT, SNP.fraction = SNP.fraction, seed = seed, BINS = BINS,SNP.test=SNP.test,DPP=DPP, SNP.permutation=SNP.permutation,
                        LD.chromosome=LD.chromosome,LD.location=LD.location,LD.range=LD.range,SNP.CV=SNP.CV,SNP.robust=SNP.robust,
                        genoFormat=genoFormat,hasGenotype=hasGenotype,byFile=byFile,fullGD=fullGD,PC=PC,GI=GI,Timmer = Timmer, Memory = Memory,
                        sangwich.top=sangwich.top,sangwich.bottom=sangwich.bottom,QC=QC,GTindex=GTindex,LD=LD,file.output=file.output,cutOff=cutOff, 
                        Model.selection = Model.selection, Create.indicator = Create.indicator,
						QTN=QTN, QTN.round=QTN.round,QTN.limit=QTN.limit, QTN.update=QTN.update, QTN.method=QTN.method)  
}# end of loop on trait

if(ncol(Y>2) &file.output)
{
Timmer=gapitMain$Timmer
Memory=gapitMain$Memory

file=paste("GAPIT.", "All",".Timming.csv" ,sep = "")
write.table(Timmer, file, quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)

file=paste("GAPIT.", "All",".Memory.Stage.csv" ,sep = "")
write.table(Memory, file, quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)
}

if(ncol(Y)==2) {
h2= as.matrix(as.numeric(as.vector(gapitMain$Compression[,5]))/(as.numeric(as.vector(gapitMain$Compression[,5]))+as.numeric(as.vector(gapitMain$Compression[,6]))),length(gapitMain$Compression[,6]),1)
colnames(h2)=c("Heritability")
  print("GAPIT accomplished successfully for single trait. Results are saved. GWAS and GPS are returned!")
  return (list(QTN=gapitMain$QTN,GWAS=gapitMain$GWAS,GPS=gapitMain$GPS,Pred=gapitMain$Pred,compression=as.data.frame(cbind(gapitMain$Compression,h2)), kinship.optimum=gapitMain$kinship.optimum,kinship=gapitMain$kinship,PCA=gapitMain$PC))
}else{
  print("GAPIT accomplished successfully for multiple traits. Results are saved")
  return (list(GWAS=NULL,GPS=NULL,Pred=NULL,compression=NULL,kinship.optimum=NULL,kinship=gapitMain$KI,PCA=gapitMain$PC))
}

}# end ofdetecting null Y

}  #end of GAPIT function

