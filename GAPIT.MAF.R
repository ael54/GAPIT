`GAPIT.MAF` <-
function(MAF=NULL,P=NULL,E=NULL,trait="",threshold.output=.01){
#Object: To display probability and effect over MAF
#Input: MAF vector of MAF
#Input: P vector of P values
#Output: A table and plot
#Requirment: NA
#Authors: Zhiwu Zhang
# Start  date: April 5, 2013 
# Last update: April 5, 2013 
##############################################################################################
#print("MAF plot started")
#print(threshold.output)
#Remove NAs and under threshold

maxColor=100
index= which(P<threshold.output & !is.na(MAF))
MAF=MAF[index] 
#E=E[index] 
P=P[index] 


LP=-log10(P) 
LPC=round(LP*10,digits = 0) 
ncolors=max(LPC,na.rm=T)

#handler of extrem samll p value
CF=1
if(ncolors>maxColor){
CF=ncolors/maxColor
ncolors=maxColor
}

#print("MAF plot started 0001")
#print(length(P))
#print(ncolors)

#palette(rainbow(ncolors))
#palette(gray(seq(.9,0,len = ncolors)))
theColor=heat.colors(ncolors, alpha = 1)
#print("MAF plot started 0001b")
palette(rev(theColor))
#print("MAF plot started 0001")

pdf(paste("GAPIT.", trait,".MAF.pdf" ,sep = ""), width = 5,height=5) 
par(mar = c(5,6,5,3))


  
plot(MAF,LP,type="p",lty = 1,lwd=2,col=LPC/CF,xlab="MAF",ylab =expression(Probability~~-log[10](italic(p))),main = trait, cex.axis=1.1, cex.lab=1.3)

#for(i in 2:nc){
#lines(power[,i]~FDR, lwd=2,type="o",pch=i,col=i)
#}
#legend("bottomright", colnames(power), pch = c(1:nc), lty = c(1,2),col=c(1:nc))

dev.off()
palette("default")      # reset back to the default

}   #GAPIT.ROC ends here

