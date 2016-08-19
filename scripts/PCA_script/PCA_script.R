devtools::install_github("vqv/ggbiplot")
library(ggbiplot)

###############################
#Set your working directory
setwd("~/Desktop/Rstudio")

#Load in data from a csv
pheno.data <- read.csv('Churchill_project/Data/svenson.phenotypes.final.850.csv',header=T)

#If a phenotype has more than na.thres NAs, it will be removed
na.thres <- 100

#Which principle components to look at (list of only 2)
princomps <- c(1,2)
###############################

pheno.data.rm <- data.frame()
for(i in 1:dim(pheno.data)[2]){
  #compare length of the variable to length of variable with NAs omitted
  if((length(pheno.data[,i]) - length(na.omit(pheno.data[,i]))) > na.thres){
    pheno.data.rm <- rbind(pheno.data.rm, i)
  }
}

pheno.data <- pheno.data[,-unlist(pheno.data.rm)]
pheno.data <- na.omit(pheno.data)

sex_by_diet <- c()
for(i in 1:nrow(pheno.data)){
  sex_by_diet <- c(sex_by_diet, 
                   paste(pheno.data[i,grep("Sex", colnames(pheno.data), ignore.case = TRUE)], 
                         pheno.data[i,grep("Diet", colnames(pheno.data), ignore.case = TRUE)], 
                         sep = "_"))
}

for(i in 1:ncol(pheno.data)){
  pheno.data[,i] <- as.numeric(as.character(pheno.data[,i]))
}

rz.transform<-function(y) {
  rankY=rank(y, ties.method="average", na.last="keep")
  rzT=qnorm(rankY/(length(na.exclude(rankY))+1))  
  rzT
}

pheno.data <- pheno.data
for(i in 1:ncol(pheno.data)){
  pheno.data[,i] <- rz.transform(pheno.data[,i])
}

pheno.data <- t(na.omit(t(pheno.data)))

pc <- princomp(pheno.data)

quartz()
ggbiplot(pc, choices = princomps, circle = TRUE, ellipse = TRUE, groups = sex_by_diet)


