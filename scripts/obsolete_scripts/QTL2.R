#devtools::install_github("simecek/intermediate")
library(qtl2scan)
library(qtl2plot)
library(qtl2geno)
library(qtl2convert)
library(intermediate)

setwd("~/Desktop/Rstudio")
#load annotations, covariates, expression, snps, kinship, and genoprobs
load("Churchill_project/Data/DO478_ExprData_4eQTLMapping.Rdata")

#Load in phenotype data
oldpheno <- read.csv('Churchill_project/Data/data_restructured.csv',header=T)
newpheno <- oldpheno[match(names(probs.478[,1,1]), oldpheno$Sample),,drop=FALSE]
rownames(newpheno) <- newpheno[,3]

#microalbumin was a factor, convert to numeric
newpheno$urine.microalbumin1 <- as.numeric(as.character(newpheno$urine.microalbumin1))

#RankZ transform the data
rz.transform<-function(y) {
  rankY=rank(y, ties.method="average", na.last="keep")
  rzT=qnorm(rankY/(length(na.exclude(rankY))+1))  
  rzT
}

#rz transform the phenotype data
newpheno.rz <- newpheno
for(i in 9:ncol(newpheno)){
  newpheno.rz[,i] <- rz.transform(newpheno.rz[,i])
}

#create covariate matrix of Sex and Diet
covar <- cbind(newpheno[,"Sex"], newpheno[,"Diet"])
colnames(covar)[1:2] = c("Sex", "Diet")
covar.gen <- cbind(covar, newpheno[,5:6])
covar.gen[,3] <- as.numeric(covar.gen[,3])-1
covar.gen[,4] <- as.numeric(covar.gen[,4])-1
colnames(covar.gen)[1:4] = c("Sex", "Diet", "Gen", "parity")


if (class(probs.478) != "calc_genoprob"){
  probs=probs_doqtl_to_qtl2(probs.478, map=snps.64K,
                            chr_column = "Chr", pos_column = "cM",
                            marker_column = "SNP_ID")
}

qtl2apr = genoprob_to_alleleprob(probs)

#get kinship and sex covars
k <- calc_kinship(qtl2apr, "loco", cores=8)
sex <- covar.gen[,1,drop = FALSE]
diet <- covar[,2,drop = FALSE]

newpheno.rz <- newpheno.rz[,9:ncol(newpheno.rz)]
pp=newpheno.rz[,c("BW.10", "BW.12"), drop = FALSE]
rownames(pp)=rownames(newpheno.rz)
colnames(pp)=names(newpheno.rz)["BW.10"]

pp.expr <- expr.rna.478[,"ENSMUSG00000055850", drop = FALSE]

#time using one core
#user  system elapsed 
#29.483   6.893  37.869 
#plot(scan1(qtl2apr, pp, k, sex, cores=4))
qtl <- scan1(qtl2apr, pp, k, addcovar = covar.gen[,2:4], intcovar = covar.gen[,1,drop=FALSE] cores = 4)

peak <- bayes_int(qtl, chr = NULL)

geneRange <- data.frame()
for(i in 1:nrow(annotations.rna.478)){
  if(annotations.rna.478[i,4] < peak[,3] && annotations.rna.478[i,5] > peak[,1] && annotations.rna.478[i,3] == 10){
    insert <- cbind(annotations.rna.478[i,2], annotations.rna.478[i,4], annotations.rna.478[i,5], i)
    geneRange <- rbind(geneRange, insert)
  }
}

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



