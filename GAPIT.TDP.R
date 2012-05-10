`GAPIT.TDP` <-
function(Y,CV=NULL,SNP){
#Object: To perform GWAS and GPS with Trillion Data Points (TDP)
#Output: GWAS
#Authors: Zhiwu Zhang
# Last update: may 9, 2012
##############################################################################################

mydataPath="C:\\myGAPIT\\"
myY  <- read.table(paste(mydataPath,"mdp_traits.txt",sep=""), head = TRUE)
myGD <- read.table(paste(mydataPath,"X_3122_SNPs.txt",sep=""), head = TRUE)
data=merge(myY[,1:2], myGD, by.x = "Taxa", by.y = "Taxa")
m=ncol(data)-2
n=nrow(data)-2
xnam <- paste("x", 1:m, sep="")
colnames(data)=c("Taxa","y",xnam)

q=floor(sqrt(n))

xnam <- paste("x", 1:m, sep="")
(fmla <- as.formula(paste("y ~ ", paste(xnam, collapse= "+"))))
myLM=lm(fmla)


names(myY)
attach(myY)

myLM=lm(data[,2]~data[,3]:data[,54])
myLM3=lm(data[,2]~data[,3])
myLM54=lm(data[,2]~data[,54])


return (list(GWAS=GWAS,Timmer=Timmer,Memory=Memory))

}#The function GAPIT.TDP ends here


