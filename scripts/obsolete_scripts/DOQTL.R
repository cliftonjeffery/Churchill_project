#code to load in DOQTL
# source("http://bioconductor.org/biocLite.R")
# biocLite("DOQTL")
library(DOQTL)
setwd("~/Desktop/Rstudio")
load("Churchill_project/Data/DO478_ExprData_4eQTLMapping.Rdata")

#Calculate Kinship
K = K.LOCO.478

#get pheno data and reduce it to only the mice in the dataset
oldpheno <- read.csv('Churchill_project/Data/data_restructured.csv',header=T)
newpheno <- oldpheno[match(names(probs.478[,1,1]), oldpheno$Sample),]

#microalbumin was a factor, convert to numeric
newpheno$urine.microalbumin1 <- as.numeric(as.character(newpheno$urine.microalbumin1))

#create covariate matrix of Sex and Diet
covar = model.matrix(~Sex, data = newpheno.rz)[,2]
covar <- cbind(covar, model.matrix(~Diet, data = newpheno.rz)[,2])
colnames(covar)[1:2] = c("Sex", "Diet")
covar.sex <- covar[,-2]
covar.sex <- as.data.frame(covar.sex)
colnames(covar.sex) <- "Sex"

#DOQTL analyzes whether a chromosome is a sex chromosome based on whether the name is numeric,
#change the X chromosome to numeric name s that is treated the same as the others
snps.64K.noX <- snps.64K
snps.64K.noX$Chr <- gsub("X", 20, snps.64K$Chr)

#QTL scan
rownames(newpheno.rz) <- rownames(probs.478)
rownames(covar) <- rownames(probs.478)
rownames(covar.sex) <- rownames(expr.rna.478)
qtl = scanone(pheno = newpheno, pheno.col = "TG1", probs = probs.478, K = K, 
              addcovar = covar.sex, snps = snps.64K.noX)
plot(qtl, main = "covar.sex")
qtl = scanone(pheno = newpheno, pheno.col = "TG1", probs = probs.478, K = K, 
              addcovar = covar, snps = snps.64K.noX)
plot(qtl, main = "covar")
pos <- bayesint(qtl, chr = 15)[2,3]
covar.gene <- cbind(covar, probs.478[,,pos])
covar.gen <- cbind(covar, newpheno.rz[,4], newpheno.rz[,5])
qtl = scanone(pheno = newpheno.rz, pheno.col = "BW.4", probs = probs.478, K = K, 
              addcovar = covar.gen, snps = snps.64K.noX)
plot(qtl, main = "BW.4")

#run permutations to find threshold values
perms <- scanone.perm(pheno = newpheno, pheno.col = "GLDH1", probs = probs.478,
                    addcovar = covar, snps = snps.64K, nperm = 2)

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
rownames(newpheno) <- rownames(probs.478)
rownames(covar) <- rownames(probs.478)
rownames(covar.sex) <- rownames(expr.rna.478)
covar.drop <- cbind(covar.sex, jitter(expr.rna.478[,17002], factor = 0.1))
colnames(covar.drop)[2] <- "Ikbke"
qtl = scanone(pheno = newpheno, pheno.col = "TG2", probs = probs.478, K = K, 
              addcovar = covar, snps = snps.64K.noX)

interval = bayesint(qtl, chr = 15)

#find all genes in significant range of qtl peak
geneRange <- data.frame()
for(i in 1:nrow(annotations.rna.478)){
  if(annotations.rna.478[i,4] < interval[3,3] && annotations.rna.478[i,5] > interval[1,3] && annotations.rna.478[i,3] == 15){
    insert <- cbind(annotations.rna.478[i,2], annotations.rna.478[i,4], annotations.rna.478[i,5], i)
    geneRange <- rbind(geneRange, insert)
  }
}

for (i in 1:dim(geneRange)[1]){
  covar.drop <- cbind(covar.sex, expr.rna.478[,grep(annotations.rna.478
            [grep(geneRange[i,1], annotations.rna.478[,2]),1],colnames(expr.rna.478))])
  qtl.drop <- scanone(pheno = expr.rna.478, pheno.col = 15878, probs = probs.478, K = K, 
                      addcovar = covar.drop, snps = snps.64K.noX)
  peak.value <- bayesint(qtl.drop, chr = 1:19)
  if(max(peak.value[,7]) < (max(interval[,7]) - 2)){
    print(geneRange[i,1])
  }
}

genoprobs <- array()
expr <- data.frame() #just get column 10
fileList <- as.vector(read.csv("folder.txt", header = FALSE))
for(i in 1:nrow(fileList)){
  #load in genoprobs data
  genoprobs.mouse <- read.table(paste(fileList[i], gbrs.interpolated.genoprobs.tsv), sep = "\t", header = FALSE)
  #set colnames
  colnames(genoprobs.mouse) <- c("A", "B", "C", "D", "E", "F", "G", "H")
  #bind the data to the already existing data
  genoprobs <- abind(genoprobs, genoprobs.mouse, along = 3)
  #load in expression data
  expr.mouse <- read.table(paste(fileList[i], gbrs.quantified.diploid.genes.expected_read_counts), sep = "\t", header = FALSE)
  #set expression data rownames to ensembl ids
  rownames(expr.mouse) <- expr.mouse[,1]
  #transpose the data so it's a single row and markers are in the columns
  totalexpr.mouse <- t(expr.mouse[,10])
  #bind the data together
  rbind(expr, expr.mouse)
}
#set mouse id to be y axis, strains to be x, and markers to be z
genoprobs <- aperm(genoprobs, c(3,1,2))
#identify mice in both sets by the folder name
rownames(genoprobs) <- fileList
rownames(expr) <- fileList

covar.peak <- covar
for(i in 1:nrow(geneRange)){
  covar.peak <- cbind(covar.peak, expr.rna.478[,geneRange[i,1]])
}

probs <- probs.478[,,interval[2,3]]
haplotypes <- c()
for (i in 1:nrow(probs)) {
  sort <- sort(probs[i,])
  haplotypes <- c(haplotypes, paste(names(sort)[7], names(sort)[8], sep = ""))
}

covar.hap <- cbind(covar, as.factor(haplotypes))

qtl = scanone(pheno = newpheno, pheno.col = "TG1", probs = probs.478, K = K, 
              addcovar = covar.peak, snps = snps.64K.noX)
plot(qtl, main = "peak")
interval <- bayesint(qtl, chr = 14)

for(i in 1:nrow(interaction)){
  if ((as.numeric(as.character(interaction[i,3])) > as.numeric(as.character(interaction[i,5]))) & 
      (as.numeric(as.character(interaction[i,4])) > as.numeric(as.character(interaction[i,5])))){
    print(as.character(interaction[i,2]))
  }
}








