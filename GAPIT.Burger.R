`GAPIT.Burger` <-
function(Y=NULL,CV=NULL,GK=NULL){
#Object: To calculate likelihood, variances and ratio  
#Straitegy: NA
#Output: P value
#intput: 
#Y: phenotype with columns of taxa,Y1,Y2...
#CV: covariate variables with columns of taxa,v1,v2...
#GK: Genotype data in numerical format, taxa goes to row and snp go to columns. the first column is taxa (same as GAPIT.bread)
#Authors: Zhiwu Zhang
#Last update: November 2, 2011
##############################################################################################

#print("GAPIT.Burger in progress...")

myFaSTREML=GAPIT.get.LL(pheno=matrix(Y[,-1],nrow(Y),1),geno=NULL,snp.pool=as.matrix(GK[,-1]),X0=as.matrix(cbind(matrix(1,nrow(CV),1),CV[,-1])))

REMLs=-2*myFaSTREML$LL  
delta=myFaSTREML$delta
vg=myFaSTREML$vg
ve=myFaSTREML$ve

#print("GAPIT.Burger succeed!")  
return (list(REMLs=REMLs,vg=vg,ve=ve,delta=delta))
} #end of GAPIT.Burger

