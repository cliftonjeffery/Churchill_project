peak <- bayes_int(allscan1.intsex[,"CHOL1"], chr = 2)

geneRange <- data.frame()
for(i in 1:nrow(annotations.rna.478)){
  if(annotations.rna.478[i,4] < peak[,3] && annotations.rna.478[i,5] > peak[,1] && annotations.rna.478[i,3] == 2){
    insert <- cbind(annotations.rna.478[i,2], annotations.rna.478[i,4], annotations.rna.478[i,5], i)
    geneRange <- rbind(geneRange, insert)
  }
}

qtl <- scan1(qtl2apr[,'2'], pp, k["2"], covar.gen)

drops <- data.frame()
for(i in 1:nrow(geneRange)){
  covar.expr <- cbind(covar.gen, expr.rna.478[,annotations.rna.478[grep(geneRange[i,1], 
                annotations.rna.478[,2]), 1]])
  qtl <- scan1(qtl2apr[,'2'], pp, k["2"], covar.expr[,2:5], intcovar = sex)
  insert <- cbind(geneRange[i,1], find_peaks(qtl)[1,5])
  drops <- rbind(drops, insert)
}
