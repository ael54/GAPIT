`GAPIT.PCA` <-
function(X,taxa, PC.number = min(ncol(X),nrow(X)),file.output=TRUE){
# Object: Conduct a principal component analysis, and output the prinicpal components into the workspace,
#         a text file of the principal components, and a pdf of the scree plot
# Authors: Alex Lipka and Hyun Min Kang
# Last update: May 31, 2011 
############################################################################################## 

#Conduct the PCA 
print("Callingg prcomp...")
PCA.X <- prcomp(X)

print("Creating PCA graphs...")
#Create a Scree plot 
if(file.output & PC.number>1) {
pdf("GAPIT.PCA.eigenValue.pdf", width = 12, height = 12)
par(mar = c(5,5,5,5))
screeplot(PCA.X, type="lines")
dev.off()

pdf("GAPIT.PCA.pdf", width = 12, height = 12)
par(mar = c(5,5,5,5))
maxPlot=min(as.numeric(PC.number[1]),3)

for(i in 1:(maxPlot-1)){
for(j in (i+1):(maxPlot)){
plot(PCA.X$x[,i],PCA.X$x[,j],xlab=paste("PC",i,sep=""),ylab=paste("PC",j,sep=""),pch=1,col="red",cex=1.0,cex.lab=1.5, cex.axis=1.2, lwd=2,las=1)
}
}
dev.off()
}

print("Joining taxa...")
#Extract number of PCs needed
PCs <- cbind(taxa,as.data.frame(PCA.X$x))

#Remove duplicate (This is taken care by QC)
#PCs.unique <- unique(PCs[,1])
#PCs <-PCs[match(PCs.unique, PCs[,1], nomatch = 0), ]

print("Exporting PCs...")
#Write the PCs into a text file
if(file.output) write.table(PCs, "GAPIT.PCA.csv", quote = FALSE, sep = ",", row.names = FALSE,col.names = TRUE)

#Return the PCs
return(list(PCs=PCs,EV=PCA.X$sdev^2,nPCs=NULL))
}

