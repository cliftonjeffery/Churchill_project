.libPaths("/home/jeffec/scripts/packages")
library(qtl2scan)
library(qtl2plot)
library(qtl2convert)
library(qtl2geno)
devtools::install_github("simecek/intermediate")
library(intermediate)
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

covar.gen <- cbind(sex, diet, as.numeric(gen[,1]) - 1, as.numeric(par[,1]) - 1)
covar.gen <- as.matrix(covar.gen)

load("../data/scans.Rdata")

colnames(annotations.rna.478) <-c(colnames(annotations.rna.478)[1:2], "CHR", "POS", 
                                  colnames(annotations.rna.478)[5:8]) 

add <- list()
intsex <- list()
intdiet <- list()
qtl.geno <- probs.478[,,rownames(max_scan1(allscan1.intsex[,"BW.20"], chr = 9))[1]]
for(i in 1:128){
 chr.add <- list()
 chr.intsex <- list()
 chr.intdiet <- list()
 for(j in 1:19){
   print(paste(colnames(newpheno.rz)[i], ", chr ", j))
   qtl.geno <- probs.478[,,rownames(max_scan1(allscan1.intsex[,i], chr = j))[1]]
   med.add <- mediation.scan(target = newpheno.rz[,i], mediator = expr.rna.478,
                         annotation = annotations.rna.478, covar = as.matrix(covar.gen), qtl.geno =
                         qtl.geno)
   chr.add[[j]] <- med.add
   names(chr.add)[j] <- j
   qtl.geno.int <- cbind(qtl.geno, qtl.geno*covar.gen[,1])
   med.add <- mediation.scan(target = newpheno.rz[,i], mediator = expr.rna.478,
                         annotation = annotations.rna.478, covar = as.matrix(covar.gen), qtl.geno =
                         qtl.geno.int)
   chr.intsex[[j]] <- med.add
   names(chr.intsex)[j] <- j
   qtl.geno.int <- cbind(qtl.geno, qtl.geno*covar.gen[,2])
   med.add <- mediation.scan(target = newpheno.rz[,i], mediator = expr.rna.478,
                         annotation = annotations.rna.478, covar = as.matrix(covar.gen), qtl.geno =
                         qtl.geno.int)
   chr.intdiet[[j]] <- med.add
   names(chr.intdiet)[j] <- j
 }
 add[[i]] <- chr.add
 intsex[[i]] <- chr.intsex
 intdiet[[i]] <- chr.intdiet
 names(add)[i] <- colnames(newpheno.rz)[i]
 names(intsex)[i] <- colnames(newpheno.rz)[i]
 names(intdiet)[i] <- colnames(newpheno.rz)[i]
}
med_scans <- list()
med_scans[[1]] <- add
med_scans[[2]] <- intsex
med_scans[[3]] <- intdiet
names(med_scans)[1] <- "add"
names(med_scans)[2] <- "intsex"
names(med_scans)[3] <- "intdiet"
print("done")
save(med_scans, file = "../results/med_scans.Rdata")
