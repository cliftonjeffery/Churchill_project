.libPaths("/home/jeffec/scripts/packages")
library(qtl2scan)
library(qtl2plot)
library(qtl2geno)
library(qtl2convert)

load("../data/DO478_ExprData_4eQTLMapping.Rdata")

rz.transform<-function(y) {
  rankY=rank(y, ties.method="average", na.last="keep")
  rzT=qnorm(rankY/(length(na.exclude(rankY))+1))
  rzT
}

oldpheno <- read.csv('../data/data_restructured.csv',header=T)
newpheno <- oldpheno[match(names(probs.478[,1,1]), oldpheno$Sample),]
rownames(newpheno) <- rownames(expr.rna.478)

#microalbumin was a factor, convert to numeric
newpheno$urine.microalbumin1 <- as.numeric(as.character(newpheno$urine.microalbumin1))
newpheno.rz <- newpheno
for(i in 9:ncol(newpheno)){
  newpheno.rz[,i] <- rz.transform(newpheno.rz[,i])
}

#create covariate matrix of Sex and Diet
covar = model.matrix(~Sex, data = newpheno.rz)[,2]
covar <- cbind(covar, model.matrix(~Diet, data = newpheno.rz)[,2])
colnames(covar)[1:2] = c("Sex", "Diet")
covar.sex <- covar[,-2]
covar.sex <- as.data.frame(covar.sex)
colnames(covar.sex) <- "Sex"
covar.gen <- cbind(covar, newpheno[,5:6])
covar.gen[,3] <- as.numeric(covar.gen[,3])
covar.gen[,4] <- as.numeric(covar.gen[,4])

if (class(probs.478) != "calc_genoprob"){
  probs=probs_doqtl_to_qtl2(probs.478, map=snps.64K,
                            chr_column = "Chr", pos_column = "cM",
                            marker_column = "SNP_ID")
}

qtl2apr = genoprob_to_alleleprob(probs)

#get kinship and sex covars
k <- calc_kinship(qtl2apr, "loco", cores=8)
sex <- covar.sex
diet <- data.frame(covar[,2])

newpheno.rz <- newpheno.rz[,9:ncol(newpheno.rz)]

pp.expr <- as.matrix(expr.rna.478[,"ENSMUSG00000032396"])

allcoef.add <- list()
allcoef.intsex <- list()
allcoef.intdiet <- list()
for(i in 1:ncol(newpheno.rz)){

  print(paste("phenotype ", i))

  pp=as.matrix(newpheno.rz[,i])
  rownames(pp)=rownames(newpheno.rz)
  colnames(pp)=names(newpheno.rz)[i]
  
  out.add <- list()
  out.intsex <- list()
  out.intdiet <- list()

  for(j in 1:19){

    print(paste("chromosome ", j))

    out1=scan1coef(qtl2apr[,as.character(j)], pp, k[[j]], addcovar = covar.gen[,2:4], intcovar = sex)
    out2=scan1coef(qtl2apr[,as.character(j)], pp, k[[j]], addcovar = covar.gen[,c(1,3:4)], intcovar = diet)
    out3=scan1coef(qtl2apr[,as.character(j)], pp, k[[j]], addcovar = covar.gen)
    out.add[[j]] <- out3
    names(out.add)[j] <- as.character(j)
    out.intsex[[j]] <- out1
    names(out.intsex)[j] <- j
    out.intdiet[[j]] <- out2
    names(out.intdiet)[j] <- j
  }
  
  allcoef.add[[i]] <- out.add
  names(allcoef.add)[i] <- colnames(newpheno.rz)[i]
  allcoef.intsex[[i]] <- out.intsex
  names(allcoef.intsex)[i] <- colnames(newpheno.rz)[i]
  allcoef.intdiet[[i]] <- out.intdiet
  names(allcoef.intdiet)[i]	<- colnames(newpheno.rz)[i]

}

save(allcoef.add, allcoef.intsex, allcoef.intdiet, file = "../results/allcoef.Rdata")
