#code to load in DOQTL
# source("http://bioconductor.org/biocLite.R")
# biocLite("DOQTL")
library(DOQTL)
setwd("~/Desktop/Rstudio")
load("Churchill_project/Data/DO478_ExprData_4eQTLMapping.Rdata")

#create contribution graph for a sample mouse
image(1:64000, 1:ncol(probs.478), t(probs.478[1,8:1,1:64000]), breaks = 0:100/100,
      col = grey(99:0/100), axes = F, xlab = "Markers", ylab = "Founders",
      main = "Founder Allele Contributions for Sample 1")
abline(h = 0:8 + 0.5, col = "grey70")
usr = par("usr")
rect(usr[1], usr[3], usr[2], usr[4])
axis(side = 1, at = 0:5 * 100, labels = 0:5 * 100)
axis(side = 2, at = 1:8, labels = LETTERS[8:1], las = 1, tick = F)

#Calculate Kinship
K = kinship.probs(probs.478, snps = snps.64K, bychr = TRUE)

#create heat map for kinship
image(1:nrow(K[[1]]), 1:ncol(K[[1]]), K[[1]][,ncol(K[[1]]):1], xlab = "Samples", 
      ylab = "Samples", yaxt = "n", main = "Kinship between samples", 
      breaks = 0:100/100, col = heat.colors(length(0:100) - 1))
axis(side = 2, at = 30 * 0:16, labels = 30 * 16:0, las = 1)

#get pheno data and reduce it to only the mice in the dataset
oldpheno <- read.csv('Churchill_project/Data/data_restructured.csv',header=T)
newpheno <- oldpheno[match(names(probs.478[,1,1]), oldpheno$Sample),]

#microalbumin was a factor, convert to numeric
newpheno$urine.microalbumin1 <- as.numeric(as.character(newpheno$urine.microalbumin1))

#create covariate matrix of Sex and Diet
covar = model.matrix(~Sex, data = newpheno)[,2]
covar <- cbind(covar, model.matrix(~Diet, data = newpheno)[,2])
colnames(covar)[1:2] = c("Sex", "Diet")
covar.sex <- covar[,-2]
covar.sex <- as.data.frame(covar.sex)
colnames(covar.sex) <- "Sex"

#DOQTL analyzes whether a chromosome is a sex chromosome based on whether the name is numeric,
#change the X chromosome to numeric name s that is treated the same as the others
snps.64K.noX <- snps.64K
snps.64K.noX$Chr <- gsub("X", 20, snps.64K$Chr)

#QTL scan
qtl = scanone(pheno = newpheno, pheno.col = "GLDH1", probs = probs.478, K = K, 
              addcovar = covar, snps = snps.64K.noX)

#run permutations to find threshold values
perms <- scanone.perm(pheno = newpheno, pheno.col = "GLDH1", probs = probs.478,
                      addcovar = covar, snps = snps.64K, nperm = 100)

#find significant thresholds
thr = get.sig.thr(perms[,,1], alpha = c(0.05, 0.1, 0.63), Xchr = FALSE)

#plot the qtl with the thresholds
plot(qtl, sig.thr = thr, sig.col = c("red", "orange", "goldenrod"), main = "GLDH1")

#dataset of all phenotypes with most significant p-values with the interaction of sex and diet
interactions <- read.csv("expr_interactions.csv")

#find indices of genes in actual data
ls <- data.frame()
for (i in 1:nrow(interactions)) {
  index <- grep(interactions[i,2], annotations.rna.478[,2])
  ls <- rbind(ls, index)
}

#find just the expression values for those genes
expr.reduced <- expr.rna.478[,ls[,1]]

#Correlate each to sex and diet
exprCorr <- data.frame()
for(i in 1:ncol(expr.reduced)){
  c <- cor(cbind(expr.reduced[,i], newpheno$Sex), use = "pairwise.complete.obs")
  c2 <- cor(cbind(expr.reduced[,i], newpheno$Diet), use = "pairwise.complete.obs")
  print(i)
  if(c[2,1] > .1 | c2[2,1] > .1){
    insert <- cbind(i, colnames(expr.reduced)[i], c[2,1], c2[2,1])
    exprCorr <- rbind(exprCorr, insert)
  }
}

#Correlate each to each other
exprCorr <- data.frame()
for(i in 1:ncol(expr.reduced)){
  for (j in 1:ncol(expr.reduced)) {
    if (i < j){
      c <- cor(cbind(expr.reduced[,i], expr.reduced[,j]), use = "pairwise.complete.obs")
      if(c[2,1] > .9 && c[2,1] != 1){
        insert <- cbind(i, j, colnames(expr.reduced)[i], colnames(expr.reduced)[j], c[2,1])
        exprCorr <- rbind(exprCorr, insert)
      }
    }
  }
}

#Run qtl scan
qtl = scanone(pheno = expr.rna.478, pheno.col = 701, probs = probs.478, K = K, 
              addcovar = covar.sex, snps = snps.64K.noX)

interval = bayesint(qtl, chr = 1)

#find all genes in significant range of qtl peak
geneRange <- data.frame()
for(i in 1:nrow(annotations.rna.478)){
  if(annotations.rna.478[i,4] < interval[3,3] && annotations.rna.478[i,5] > interval[1,3] && annotations.rna.478[i,3] == 1){
    insert <- cbind(annotations.rna.478[i,2], annotations.rna.478[i,4], annotations.rna.478[i,5], i)
    geneRange <- rbind(geneRange, insert)
  }
}

