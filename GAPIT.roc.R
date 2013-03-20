`GAPIT.ROC` <-
function(t=NULL,se=NULL,Vp=1,trait=""){
#Object: To make table and plot for ROC (power vs FDR)
#Input: t and se are the vectors of t value and their standard error
#Input: Vp is phenotypic variance and trait is name of the phenotype
#Output: A table and plot
#Requirment: error df is same for all SMP or large
#Authors: Zhiwu Zhang
# Last update: Feb 11, 2013 
##############################################################################################

#test
#n=1000
#trait="test"
#t=rnorm(n)
#se=sqrt(abs(rnorm(n))  )
#Vp=10

#Remove NAs
index=is.na(t)
t=t[!index]
se=se[!index]
#print(head(cbind(t,se)))
#Configration
FDR=c(0,.01,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1)
coefficient=c(0,0.01,.02,.05,.1,.2,.3)

#Power holder
nf=length(FDR)
nc=length(coefficient)
power=matrix(NA,nf,nc)

#Handler of matrix format
if(length(t)==1)t=t[,1]
if(length(se)==1)se=se[,1]

n=length(t)

#Discard negative
t=abs(t)

#sort t and se
position=order(t,decreasing = TRUE)
t=t[position]
se=se[position]
EFFECT=coefficient*sqrt(Vp)
newbit=matrix(1/se,n,1)%*%EFFECT   #n by nc matrix
tnew=newbit+t  #n by nc matrix

for (i in 1:nf){
fdr=FDR[i]
cutpoint=floor(n*fdr)
cutoff=t[cutpoint]


for (j in 1:nc){
effect= EFFECT[j]
singnificant=tnew[,j]>cutoff
count=length(t[singnificant])
power[i,j]=count/n

} #end of for on fdr
} #end of for on effect

#output
rownames(power)=FDR
colnames(power)=paste("QTN=",coefficient,sep="")
write.table(power,file=paste("GAPIT.",trait,".ROC.csv",sep=""),quote = TRUE, sep = ",", row.names = TRUE,col.names = NA)

palette(c("black","red","blue","brown", "orange","cyan", "green",rainbow(nc)))

pdf(paste("GAPIT.", trait,".ROC.pdf" ,sep = ""), width = 8,height=8) 
par(mar = c(5,5,5,5))
  
plot(FDR,power[,1],type="o",lwd=2,col=1,ylab="Power")
for(i in 2:nc){
lines(power[,i]~FDR, lwd=2,type="o",pch=i,col=i)
}
legend("bottomright", colnames(power), pch = c(1:nc), lty = c(1,2),col=c(1:nc))

dev.off()

}   #GAPIT.ROC ends here

