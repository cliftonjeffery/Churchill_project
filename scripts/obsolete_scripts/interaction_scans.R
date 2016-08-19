setwd("~/Desktop/Rstudio")
library(qtl2scan)
library(qtl2plot)
library(qtl2geno)
library(qtl2convert)
load("Churchill_project/Data/DO478_ExprData_4eQTLMapping.Rdata")

print("function")
rz.transform<-function(y) {
  rankY=rank(y, ties.method="average", na.last="keep")
  rzT=qnorm(rankY/(length(na.exclude(rankY))+1))
  rzT
}

print("read in data")

oldpheno <- read.csv('Churchill_project/Data/data_restructured.csv',header=T)
newpheno <- oldpheno[match(names(probs.478[,1,1]), oldpheno$Sample),]
rownames(newpheno) <- rownames(expr.rna.478)

#microalbumin was a factor, convert to numeric
newpheno$urine.microalbumin1 <- as.numeric(as.character(newpheno$urine.microalbumin1))

print("rz transform")

newpheno.rz <- newpheno
for(i in 9:ncol(newpheno)){
  newpheno.rz[,i] <- rz.transform(newpheno.rz[,i])
}

print("covar")

#create covariate matrix of Sex and Diet
covar = model.matrix(~Sex, data = newpheno.rz)[,2]
covar <- cbind(covar, model.matrix(~Diet, data = newpheno.rz)[,2])
colnames(covar)[1:2] = c("Sex", "Diet")
covar.sex <- covar[,-2]
covar.sex <- as.data.frame(covar.sex)
colnames(covar.sex) <- "Sex"
sex <- covar.sex
diet <- as.data.frame(covar[,2])
gen <- as.data.frame(newpheno[,5])
par <- as.data.frame(newpheno[,6])
rownames(sex) <- rownames(newpheno.rz)
rownames(diet) <- rownames(newpheno.rz)
rownames(gen) <- rownames(newpheno.rz)
rownames(par) <- rownames(newpheno.rz)

print("probs")

if (class(probs.478) != "calc_genoprob"){
  probs=probs_doqtl_to_qtl2(probs.478, map=snps.64K,
                            chr_column = "Chr", pos_column = "cM",
                            marker_column = "SNP_ID")
}

qtl2apr = genoprob_to_alleleprob(probs)


print("kinship")

#get kinship and sex covars
k <- calc_kinship(qtl2apr, "loco", cores=8)
sex <- covar.sex
diet <- data.frame(covar[,2])

print("calc pp")

newpheno.rz <- newpheno.rz[,9:ncol(newpheno.rz)]
pp=as.matrix(newpheno.rz)
rownames(pp)=rownames(newpheno.rz)
colnames(pp)=names(newpheno.rz)

covar <- cbind(sex, diet, as.numeric(gen[,1]) - 1, as.numeric(par[,1]) - 1)
covar <- as.matrix(covar)

print("add")

allscan1.add=scan1(qtl2apr, pp, k, addcovar = covar, cores=8)

print("intsex")

allscan1.intsex=scan1(qtl2apr, pp, k, covar[,c(1, 3:4)], intcovar = covar[,2], cores=3)

print("intdiet")

allscan1.intdiet=scan1(qtl2apr, pp, k, (sex + gen + par), intcovar = diet, cores=3)

print("save")

save(allscan1.add, allscan1.intdiet, allscan1.intsex, file = "../results/scans.RData")

print("dif")

allscan1.difsex <- allscan1.intsex
allscan1.difsex$lod <- allscan1.intsex$lod - allscan1.add$lod

allscan1.difdiet <- allscan1.intdiet
allscan1.difdiet$lod <- allscan1.intdiet$lod - allscan1.add$lod
save(allscan1.difsex, allscan1.difdiet, file = "../results/diffs.RData")