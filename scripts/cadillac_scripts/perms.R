.libPaths("/home/jeffec/scripts/packages")
library(qtl2scan)
library(qtl2plot)
library(qtl2convert)
library(qtl2geno)
#install.packages("StatRank", repos="http://cran.rstudio.com/")
library(StatRank)

load("../data/DO478_ExprData_4eQTLMapping.Rdata")

print("function")
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
pp=as.matrix(newpheno.rz[,"TTM2"])
rownames(pp)=rownames(newpheno.rz)
colnames(pp)=names(newpheno.rz)["TTM2"]

pp[,1] <- scramble(pp[,1])

covar <- cbind(sex, diet, as.numeric(gen[,1]) - 1, as.numeric(par[,1]) - 1)
covar <- as.matrix(covar)

perm_peaks <- data.frame()
for(i in 1:1000){
    pp[,1] <- scramble(pp[,1])
    print(i)
    allscan1.int=scan1(qtl2apr, pp, k, addcovar = covar[,2:4], intcovar = as.matrix(covar[,1]), cores = 28)
    allscan1.add <- scan1(qtl2apr, pp, k, addcovar = covar, cores = 28)
    allscan1.dif <- allscan1.add
    allscan1.dif$lod <- allscan1.int$lod - allscan1.add$lod
    perm_peaks <- rbind(perm_peaks, max(find_peaks(allscan1.dif)[,5]))
}

print("done")
write.csv(perm_peaks, "../results/perm_peaks_difdiet.csv")
