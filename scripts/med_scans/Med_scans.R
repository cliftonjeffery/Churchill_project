library(intermediate)
library(qtl2scan)
#################################
mediator <- expr.rna.478

pheno.data <- read.csv('Churchill_project/Data/svenson.phenotypes.final.850.csv',header=T)

phenotype <- "BW.20"

chr <- 9

annotations <- annotations.rna.478

map <- snps.64K

probs <- probs.478

covar <- covar.gen

addcovar <- c("Diet", "Gen", "parity")

intcovar <- "Sex"
#################################

rz.transform<-function(y) {
  rankY=rank(y, ties.method="average", na.last="keep")
  rzT=qnorm(rankY/(length(na.exclude(rankY))+1))  
  rzT
}

phenotype <- as.data.frame(pheno.data[,phenotype])

if("Sample" %in% colnames(pheno.data)){
  rownames(phenotype) <- pheno.data$Sample
}

phenotype <- as.data.frame(phenotype[which(rownames(phenotype) %in% rownames(mediator)),,drop = FALSE])

phenotype[,1] <- rz.transform(phenotype[,1])

if (class(probs.478) != "calc_genoprob"){
  probs=probs_doqtl_to_qtl2(probs, map=map,
                            chr_column = "Chr", pos_column = "cM",
                            marker_column = "SNP_ID")
}

qtl2apr = genoprob_to_alleleprob(probs)

k <- calc_kinship(qtl2apr, "loco", cores=8)

qtl <- scan1(qtl2apr, phenotype, k, addcovar = covar[,addcovar], intcovar = covar[,intcovar, drop = FALSE])

qtl.geno <- probs.478[,,rownames(max_scan1(qtl, chr = chr))[1]]

if(!is.null(intcovar)){
  qtl.geno <- cbind(qtl.geno, qtl.geno*covar[,intcovar])
}

phenotype[,1] <- as.numeric(phenotype[,1])

covar <- covar[,c(addcovar, intcovar)]

for(i in 1:ncol(covar)){
  covar[,i] <- as.numeric(covar[,i])
}

index <- grep("chr", colnames(annotations), ignore.case = TRUE)
colnames(annotations)[index] <- "CHR"
index <- grep("Start.Mbp", colnames(annotations), ignore.case = TRUE)
colnames(annotations)[index] <- "POS"
index <- grep("pos", colnames(annotations), ignore.case = TRUE)
colnames(annotations)[index] <- "POS"

med.add <- mediation.scan(target = phenotype[,1], mediator = mediator,
                          annotation = annotations, 
                          covar = as.matrix(covar), 
                          qtl.geno = qtl.geno)

kplot(med.add)
plot(med.add)
