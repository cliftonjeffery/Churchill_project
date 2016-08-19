library(DOQTL)
setwd("~/Desktop/Rstudio")
load("Churchill_project/Data/DO478_ExprData_4eQTLMapping.Rdata")
K = kinship.probs(probs.478, snps = snps.64K, bychr = TRUE)
oldpheno <- read.csv('Churchill_project/Data/data_restructured.csv',header=T)
newpheno <- oldpheno[match(names(probs.478[,1,1]), oldpheno$Sample),]
newpheno$urine.microalbumin1 <- as.numeric(as.character(newpheno$urine.microalbumin1))
covar = model.matrix(~Sex, data = newpheno)[,2]
covar <- cbind(covar, model.matrix(~Diet, data = newpheno)[,2])
colnames(covar)[1:2] = c("Sex", "Diet")
covar.sex <- covar[,-2]
covar.sex <- as.data.frame(covar.sex)
snps.64K.noX <- snps.64K
snps.64K.noX$Chr <- gsub("X", 20, snps.64K$Chr)

pdf("pheno_qtl.pdf")
bayes <- data.frame()
for (i in 9:ncol(newpheno)){
  qtl = scanone(pheno = newpheno, pheno.col = i, probs = probs.478, K = K, 
                addcovar = covar.sex, snps = snps.64K.noX)
  perms = scanone.perm(pheno = newpheno, pheno.col = i, probs = probs.478,
                       addcovar = covar.sex, snps = snps.64K.noX, nperm = 100)
  thr = get.sig.thr(perms[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = TRUE)
  plot(qtl, sig.thr = thr, sig.col = c("red", "orange", "goldenrod"), main = colnames(newpheno)[i])
  for(j in 1:20){
    if (max(bayesint(qtl, chr = j)[,7]) > thr[,2]){
      insert <- cbind(colnames(newpheno)[i], paste("chr ", as.character(j)), max(bayesint(qtl, chr = j)[,7]))
      bayes <- rbind(bayes, insert)
    }
  }
}
graphics.off()

pdf("expr_qtl.pdf")
bayes_expr <- data.frame()
for (i in 1:ncol(expr.rna.478)){
  qtl = scanone(pheno = expr.rna.478, pheno.col = i, probs = probs.478, K = K, 
                addcovar = covar.sex, snps = snps.64K.noX)
  perms = scanone.perm(pheno = newpheno, pheno.col = i, probs = probs.478,
                       addcovar = covar.sex, snps = snps.64K.noX, nperm = 100)
  thr = get.sig.thr(perms[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = TRUE)
  plot(qtl, sig.thr = thr, sig.col = c("red", "orange", "goldenrod"), main = colnames(newpheno)[i])
  for(j in 1:20){
    if (max(bayesint(qtl, chr = j)[,7]) > thr[,2]){
      insert <- cbind(colnames(newpheno)[i], paste("chr ", as.character(j)), max(bayesint(qtl, chr = j)[,7]))
      bayes_expr <- rbind(bayes_expr, insert)
    }
  }
}
graphics.off()