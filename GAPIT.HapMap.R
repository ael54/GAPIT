`GAPIT.HapMap` <-
function(G,SNP.effect="Add",SNP.impute="Middle",heading=TRUE){
#Object: To convert character SNP genotpe to numerical
#Output: Coresponding numerical value
#Authors: Feng Tian and Zhiwu Zhang
# Last update: May 30, 2011 
##############################################################################################

print(paste("Converting hampmap format to numerical under model of ", SNP.impute,sep=""))
#gc()
#GAPIT.Memory.Object(name.of.trait="HapMap.Start")

#GT=data.frame(G[1,-(1:11)])
if(heading){
GT= t(G[1,-(1:11)])
GI= G[-1,c(1,3,4)]
}else{
GT=NULL
GI= G[,c(1,3,4)]
}


#Set column names
if(heading)colnames(GT)="taxa"
colnames(GI)=c("SNP","Chromosome","Position")

#Initial GD
GD=NULL
bit=nchar(as.character(G[2,12])) #to determine number of bits of genotype
#print(paste("Number of bits for genotype: ", bit))

print("Perform numericalization")
if(heading){
  GD= apply(G[-1,-(1:11)],1,function(one) GAPIT.Numericalization(one,bit=bit,effect=SNP.effect,impute=SNP.impute))
}else{
  GD= apply(G[  ,-(1:11)],1,function(one) GAPIT.Numericalization(one,bit=bit,effect=SNP.effect,impute=SNP.impute))
}

#set GT and GI to NULL in case of null GD
if(is.null(GD)){
  GT=NULL
  GI=NULL
}

print(paste("Succesfuly finished converting hampmap which has bits of ", bit,sep=""))
return(list(GT=GT,GD=GD,GI=GI))
}#end of GAPIT.HapMap function

