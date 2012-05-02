`GAPIT.Create.Indicator` <-
function(xs ){
#Object: To esimate variance component by using EMMA algorithm and perform GWAS with P3D/EMMAx
#Output: ps, REMLs, stats, dfs, vgs, ves, BLUP,  BLUP_Plus_Mean, PEV
#Authors: Alex Lipka
# Last update: April 30, 2012
##############################################################################################

unique.sorted <- unique(xs)[order(unique(xs))]
x.min <- which(unique.sorted == min(unique(xs)))
 

x.ind <- NULL


for(i in unique.sorted[-x.min]){
 x.col <- rep(NA, length(xs))
 x.col[which(xs==i)] <- 1
 x.col[which(xs!=i)] <- 0
 x.col[which(is.na(xs))] <- NaN
 x.ind <- cbind(x.ind,x.col)                         
}



return(x.ind)

#print("GAPIT.EMMAxP3D accomplished successfully!")
}#end of GAPIT.EMMAxP3D function





