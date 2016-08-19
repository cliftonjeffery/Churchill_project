#devtools::install_github("simecek/intermediate")
library(qtl2scan)
library(qtl2plot)
library(qtl2geno)
library(qtl2convert)
library(intermediate)

setwd("~/Desktop/Rstudio")
#load annotations, covariates, expression, snps, kinship, and genoprobs
load("Churchill_project/Data/DO478_ExprData_4eQTLMapping.Rdata")

#RankZ transform the data
rz.transform<-function(y) {
  rankY=rank(y, ties.method="average", na.last="keep")
  rzT=qnorm(rankY/(length(na.exclude(rankY))+1))  
  rzT
}

#Load in phenotype data
oldpheno <- read.csv('Churchill_project/Data/data_restructured.csv',header=T)

#Remove all mice that we didn't have expression data for
newpheno <- oldpheno[match(names(probs.478[,1,1]), oldpheno$Sample),]
rownames(newpheno) <- rownames(expr.rna.478)

#microalbumin was a factor, convert to numeric
newpheno$urine.microalbumin1 <- as.numeric(as.character(newpheno$urine.microalbumin1))

#rz transform the phenotype data
newpheno.rz <- newpheno
for(i in 9:ncol(newpheno)){
  newpheno.rz[,i] <- rz.transform(newpheno.rz[,i])
}

#create covariate matrix of Sex and Diet
covar = model.matrix(~Sex, data = newpheno.rz) [,2]
covar <- cbind(covar, model.matrix(~Diet, data = newpheno.rz)[,2])
colnames(covar)[1:2] = c("Sex", "Diet")
covar.sex <- covar[,-2]
covar.sex <- as.data.frame(covar.sex)
colnames(covar.sex) <- "Sex"
covar.gen <- cbind(covar, newpheno[,5:6])
covar.gen[,3] <- as.numeric(covar.gen[,3])-1
covar.gen[,4] <- as.numeric(covar.gen[,4])-1


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
pp=as.matrix(newpheno.rz[,"BW.10"])
rownames(pp)=rownames(newpheno.rz)
colnames(pp)=names(newpheno.rz)["BW.10"]

pp.expr <- as.matrix(expr.rna.478[,"ENSMUSG00000055850"])

#time using one core
#user  system elapsed 
#29.483   6.893  37.869 
#plot(scan1(qtl2apr, pp, k, sex, cores=4))
allscan1.add=scan1(qtl2apr, pp, k, covar.gen, cores=8)
allscan1.intsex=scan1(qtl2apr, pp, k, diet, intcovar = sex, cores=3)
allscan1.intdiet=scan1(qtl2apr, pp, k, sex, intcovar = diet, cores=3)

allscan1.difsex <- allscan1.intsex
allscan1.difsex$lod <- allscan1.intsex$lod - allscan1.add$lod

allscan1.difdiet <- allscan1.intdiet
allscan1.difdiet$lod <- allscan1.intdiet$lod - allscan1.add$lod

save(allscan1.add, allscan1.difsex, allscan1.difdiet, allscan1.intdiet, allscan1.intsex, file = "Churchill_project/results/scans.RData")

peak <- bayes_int(allscan1.difdiet[,"BW.4"], chr = 7)

geneRange <- data.frame()
for(i in 1:nrow(annotations.rna.478)){
  if(annotations.rna.478[i,4] < peak[,3] && annotations.rna.478[i,5] > peak[,1] && annotations.rna.478[i,3] == 10){
    insert <- cbind(annotations.rna.478[i,2], annotations.rna.478[i,4], annotations.rna.478[i,5], i)
    geneRange <- rbind(geneRange, insert)
  }
}

colnames(annotations.rna.478) <-c(colnames(annotations.rna.478)[1:2], "CHR", "POS", 
                                  colnames(annotations.rna.478)[5:8]) 

med.scans <- list()
for(i in 1:128){
  med.add <- mediation.scan(target = newpheno.rz[,i], mediator = expr.rna.478, 
                        annotation = annotations.rna.478, covar = as.matrix(covar.gen), qtl.geno = 
                        probs.478[,,rownames(max_scan1(allscan1.add[,i]))[1]])
  med.scans[[i]] <- med.add
  names(med.scans)[i] <- colnames(newpheno.rz)[i]
}

qtl <- scan1(qtl2apr, pp.expr, k, covar.gen, cores = 4)

med.expr <- mediation.scan(target = expr.rna.478[,"ENSMUSG00000038370"], mediator = expr.rna.478, 
                          annotation = annotations.rna.478, covar = cbind(covar.gen[,1], 
                          covar.gen[,3], covar.gen[,3], covar.gen[,4]), qtl.geno = 
                          probs.478[,,rownames(qtl$lod)[grep(max(qtl$lod), qtl$lod)]])


for (i in 1:dim(geneRange)[1]){
  covar.drop <- cbind(covar.gen, expr.rna.478[,annotations.rna.478
                                                    [grep(geneRange[i,1], annotations.rna.478[,2]),1]])
  qtl.drop <- scan1(qtl2apr, pp, k, 
                      addcovar = covar.drop[,c(1, 3:5)], intcovar = diet, cores = 7)
  peak.value <- find_peaks(qtl.drop)
  if(peak.value[17,5] < 10.661857){
    print(geneRange[i,1])
  }
}

Dis3l <- expr.rna.478[,"ENSMUSG00000032396"]
covar.Dis3l <- cbind(covar.gen, Dis3l)
covar.Dis3l <- as.matrix(covar.Dis3l)

scan_on_Dis3l <- scan1(qtl2apr, pp, k, covar.Dis3l, cores=3)

for(i in 1:128){
  if(med.scans[[i]]$LOD[11277] < (mean(med.scans[[i]]$LOD) * 0.8)){
    print(names(med.scans)[i])
  }
}



