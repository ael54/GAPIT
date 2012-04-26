`GAPIT.Specify` <-
function(GI=NULL,GP=NULL,bin.size=10000000,inclosure.size=NULL,MaxBP=1e10){
#Object: To get indicator (TURE or FALSE) for GI based on GP 
#Straitegy
#       1.set bins for all snps in GP
#       2.keep the snp with smallest P value in each bin, record SNP ID
#       3.Search GI for SNP with SNP ID from above
#       4.return the position for SNP selected
#Input:
#GI: Data frame with three columns (SNP name, chr and base position)
#GP: Data frame with seven columns (SNP name, chr and base position, P, MAF,N,effect)
#Output: 
#theIndex: a vector indicating if the SNPs in GI belong to QTN or not)
#Authors: Zhiwu Zhang
#Last update: September 24, 2011
##############################################################################################

#print("Specification in process...")
if(is.null(GP))return (list(index=NULL,BP=NULL))

#set inclosure bin in GP

#Create SNP ID: position+CHR*MaxBP
ID.GP=as.numeric(as.vector(GP[,3]))+as.numeric(as.vector(GP[,2]))*MaxBP

#Creat bin ID
bin.GP=floor(ID.GP/bin.size )

#Create a table with bin ID, SNP ID and p value (set 2nd and 3rd NA temporately)
binP=as.matrix(cbind(bin.GP,NA,NA,ID.GP,as.numeric(as.vector(GP[,4])))  )
n=nrow(binP)

#Sort the table by p value and then bin ID (e.g. sort p within bin ID)
binP=binP[order(as.numeric(as.vector(binP[,5]))),]  #sort on P alue
binP=binP[order(as.numeric(as.vector(binP[,1]))),]  #sort on bin

#set indicator (use 2nd 3rd columns)
binP[2:n,2]=binP[1:(n-1),1]
binP[1,2]=0 #set the first
binP[,3]= binP[,1]-binP[,2]


#Se representives of bins
ID.GP=binP[binP[,3]>0,]


#Choose the most influencial bins as estimated QTNs
ID.GP=ID.GP[order(as.numeric(as.vector(ID.GP[,5]))),]  #sort on P alue
ID.GP=ID.GP[!is.na(ID.GP[,4]),4] #must have chr and bp information, keep SNP ID only

if(!is.null(inclosure.size)) ID.GP=ID.GP[1:inclosure.size] #keep the top ones selected

#create index in GI
theIndex=NULL
if(!is.null(GI)){
ID.GI=as.numeric(as.vector(GI[,3]))+as.numeric(as.vector(GI[,2]))*MaxBP
theIndex=ID.GI %in% ID.GP
}
#print("Specification in process done")
return (list(index=theIndex,CB=ID.GP))
} #end of GAPIT.Specify

