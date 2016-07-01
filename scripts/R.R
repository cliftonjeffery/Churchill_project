library(ggplot2)
library(corrplot)
library(plyr)
library(dplyr)
library(pvclust)
library(devtools)
setwd("~/Desktop/Rstudio")

# my.data <- read.csv('Churchill_project/Data/svenson.phenotypes.final.850.csv',header=T)
# 
# #my.data <- my.data[,-c(17:44)]
# 
# #Convert Sex and Diet to numeric variables
# my.data$Sex <- gsub("F", 1, my.data$Sex)
# my.data$Sex <- gsub("M", 0, my.data$Sex)
# my.data$Sex <- as.integer(my.data$Sex)
# 
# my.data$Diet <- gsub("hf", 1, my.data$Diet)
# my.data$Diet <- gsub("chow", 0, my.data$Diet)
# my.data$Diet <- as.integer(my.data$Diet)
# 
# index <- grep("G5", my.data$Gen)
# my.data$Gen[index] <- "G5L1"
# #my.data$Gen <- gsub("G5", "G5L1", my.data$Gen)
# ls <- list(1:846)
# my.data <- cbind(my.data, ls)
# for(i in 1:dim(my.data)[1]){
#   x <-  toString(substring(my.data$Gen[i], unlist(gregexpr("L", my.data$Gen[i])), nchar(my.data$Gen[i])))
#   my.data[i,151] <- x
# }
# colnames(my.data) <- c(colnames(my.data[,1:150]), "Parity")
# my.data <- my.data[c(1:4, 151, 5:150)]
# for(i in 1:dim(my.data)[1]){
#   x <- substring(my.data$Gen[i], 1, unlist(gregexpr("L", my.data$Gen[i]))-1)
#   my.data[i,4] <- x
# }
# 
# my.data$Sex <- as.factor(my.data$Sex)
# 
# ls <- subset(naList, as.numeric(as.character(naList[,2])) > 300)
# index=as.numeric(match(as.character(ls[,1]),as.character(names(my.data))))
# my.data <- my.data[,-index]

my.data <- read.csv('Churchill_project/Data/data_restructured.csv',header=T)

my.data[,c(4,7,9:136)] <- as.data.frame(matrix(as.numeric(as.matrix(my.data[,c(4,7,9:136)])), nrow=nrow(my.data[,c(4,7,9:136)])))
#colnames(datmat) <- names(my.data[,c(4,7,9:136)])

my.data$urine.microalbumin1 <- as.numeric(my.data$urine.microalbumin1)
my.data$urine.microalbumin2 <- as.numeric(my.data$urine.microalbumin2)

#subset the data into datasets of all males, all hf, etc.
my.data.f <- subset(my.data, my.data$Sex == 1)
my.data.m <- subset(my.data, my.data$Sex == 0)

my.data.chow <- subset(my.data, my.data$Diet == 0)
my.data.hf <- subset(my.data, my.data$Diet == 1)

#Find all phenotypes that have a .05 p-value with Sex, Diet, or something else

#Significant####################

#Variable Values
pvalue_threshold <- .05
tData <- my.data
################################

significant <- data.frame()
for (i in (grep("Coat.Color", colnames(data))+1):dim(data)[2]) {
  #Finding linear models for all the variables influence on a phenotype and various subsets of the variables
  test1 <- lm(tData[,i]~tData$Sex+tData$Diet+tData$Gen + tData$parity)
  test2 <- lm(tData[,i]~tData$Sex+tData$Gen + tData$parity)
  test3 <- lm(tData[,i]~tData$Diet+tData$Gen + tData$parity)
  test4 <- lm(tData[,i]~tData$Gen + + tData$parity)
  test5 <- lm(tData[,i]~tData$Sex*tData$Diet + tData$Sex + tData$Diet + tData$Gen + tData$parity)
  test6 <- lm(tData[,i]~tData$Sex+tData$Diet+tData$Gen)
  #tests for na
  if(!is.na(anova(test1, test2)[2,6]) & !is.na(anova(test1, test3)[2,6]) & !is.na(anova(test1, test5)[2,6] < .05)){
    #tests if p-value is significant
    if(anova(test1, test2)[2,6] < pvalue_threshold | anova(test1, test3)[2,6] < pvalue_threshold | anova(test1, test5)[2,6] < pvalue_threshold | anova(test1, test6)[2,6] < pvalue_threshold){
      #create dataframe of p-values for the influence of just sex, just diet, etc. on the phenotype
      insert <- cbind(colnames(tData[i]), anova(test1, test2)[2,6], anova(test1, test3)[2,6], anova(test1, test5)[2,6], anova(test1, test6)[2,6], i, colnames(data)[i])
      significant=rbind(significant, insert)
    }
  }
}

for(i in 1:length(significant[,6])){
  if(length(grep(significant[i,6], significant.p[,6])) == 0){
    print(toString(significant[i,7]))
  }
}

#transform to Z-scores####

#Rank-z transform the data
my.data.b <- my.data
for(i in 7:dim(my.data)[1]){
  my.data[,i] <- scale(my.data[,i])
}

my.data.rz <- my.data
my.data <- my.data.b

my.data.f.rz <- subset(my.data.rz, my.data$Sex == 1)
my.data.m.rz <- subset(my.data.rz, my.data$Sex == 0)


#correlations (Obsolete)######

#Variable Values
corr_phenotype <- "Diet"
corr_threshold <- .4
##############################

#finds correlation values between diet and all phenotypes that Diet and/or Sex have a p-value less
#than .05 with (nonfunctional)

#Obsolete

correlations <- data.frame()
index <- grep(corr_phenotype, colnames(my.data))
#loops through all phenotypes with at least one significant p-value found in the above process
for (i in 1:dim(significant)[1]){
  #finds correlation between input variable and each variable in significant
  #index values for my.data are stored in significant for each variable in column 5
  mycor <- cor(cbind(as.numeric(my.data[,index]), as.numeric(my.data[,as.numeric(significant[i,5])])), use = "pairwise.complete.obs")
  if (!is.na(mycor[2,1])){
    #if the absolute value of the correlation is greater than the input threshold, add to dataframe
    if(as.numeric(mycor[2,1]) > corr_threshold | as.numeric(mycor[2,1]) < -1*corr_threshold ){
      insert <- cbind(toString(significant[i,1]), mycor[2,1], toString(significant[i,3]), i)
      correlations=rbind(correlations, insert)
    }
  }
}
dim(correlations)

#Find correlations that go opposite ways

#Non-functional (yet)

oppCorr <-data.frame()
for(i in (grep("Coat.Color", colnames(my.data))+1):dim(my.data)[2]){
  for(j in (grep("Coat.Color", colnames(my.data))+1):dim(my.data)[2]){
    mycor.hf <- cor(cbind(as.numeric(my.data.hf[,i]), as.numeric(my.data.hf[,j])), use = "pairwise.complete.obs")
    mycor.chow <- cor(cbind(as.numeric(my.data.chow[,i]), as.numeric(my.data.chow[,j])), use = "pairwise.complete.obs")
    if (!is.na(mycor.hf[2,1]) &! is.na(mycor.chow[2,1])){
      if(as.numeric((mycor.hf[2,1]) > 0 | as.numeric(mycor.hf[2,1]) < 0) & as.numeric(mycor.hf[2,1]) != 1){
        if(as.numeric((mycor.chow[2,1]) > 0 | as.numeric(mycor.chow[2,1]) < 0) & as.numeric(mycor.chow[2,1]) != 1){
          print(" ")
          if((mycor.hf > 0 & mycor.chow < 0) | (mycor.hf < 0 & mycor.chow > 0)) {
            print(" ")
            insert <- cbind(colnames(my.data[i]), colnames(my.data[j]), mycor.hf[2,1], mycor.chow[2,1])
            oppCorr <- rbind(oppCorr, insert)
          }
        }
      }
    }
  }
}

#allCorr####################

#Variable Values
sig_corr <- .3
tData <- my.data
############################

#Correlate everything together
allCorr <-data.frame()
#leaving out generaion, coat color, and sample
for(i in c(grep("Sex", colnames(tData)),grep("Diet", colnames(tData)),grep("ACR1", colnames(tData)):dim(tData)[2])){
  for(j in c(grep("Sex", colnames(tData)),grep("Diet", colnames(tData)),grep("ACR1", colnames(tData)):dim(tData)[2])){
    #Eliminates duplicate entries
    if (i < j){
      #find the correlation between all variables in dataset
      mycor <- cor(cbind(as.numeric(tData[,i]), as.numeric(tData[,j])), use = "pairwise.complete.obs")
      if (!is.na(mycor[2,1])){
        #If correlation coefficient above selected threshold, store value
        if(as.numeric((mycor[2,1]) > sig_corr | as.numeric(mycor[2,1]) < -sig_corr) & as.numeric(mycor[2,1]) != 1){
          insert <- cbind(colnames(tData[i]), colnames(tData[j]), mycor[2,1])
          allCorr <- rbind(allCorr, insert)
        }
      }
    }
  }
}

#betterCorr#########################

#Variable Values
phenotype <- "Diet" #Case sensitive
####################################

#Find all significant variables related to a certain phenotype
betterCorr <- data.frame()
#Check each correlation
for (i in 1:dim(allCorr)[1]){
  #if input variable is in entry and the phenotype it is correlated to has
  #a significant p-value in one of the three tests done in the first loop
  if((allCorr[i,1] == phenotype | allCorr[i,2] == phenotype) && (length(grep(allCorr[i,2], significant[,1])) != 0)){
    #Add this correlation value to a data-frame
    insert <- cbind(toString(allCorr[i,2]), toString(allCorr[i,3]), toString(significant[grep(allCorr[i,2], significant[,1]),2]))
    betterCorr <- rbind(betterCorr, insert)
  }
}

#bwCorr############################

#check if HF diet is positively correlated with Body Weight
correlations.BW <- data.frame()
for (i in c(7:16, 45:149)){
  mycor <- cor(cbind(as.numeric(my.data[,150]), as.numeric(my.data[,i])), use = "pairwise.complete.obs")
  if (!is.na(mycor[2,1])){
    if(as.numeric((mycor[2,1]) > .4 | as.numeric(mycor[2,1]) < -.4) & as.numeric(mycor[2,1]) != 1){
      insert <- cbind(colnames(my.data[i]), mycor[2,1])
      correlations.BW=rbind(correlations.BW, insert)
    }
  }
}

#sigCorr##############

#Variable Values
pval_threshold <- 0.05
######################

#Find all phenotypes correlated to each other, with both having significant p-values
#Phenoype is changed in betterCorr
sigCorr <- data.frame()
#Run through every correlation found in allCorr
for(i in 1:dim(allCorr)[1]){
  #Check if each correlated variable is in significant
  if(length(significant[grep(allCorr[i,1], significant[,1]), 1]) != 0 && length(significant[grep(allCorr[i,2], significant[,1]), 1]) != 0){
    #Check if the pvalue is NA
    if(!is.na(as.numeric(toString(significant[grep(allCorr[i,1], significant[,1]), 2])))  && !is.na(as.numeric(toString(significant[grep(allCorr[i,2], significant[,1]), 2])))){
      #check if the pvalue is less than the threshold
      if(as.numeric(toString(significant[grep(allCorr[i,1], significant[,1]), 2])) < pval_threshold && as.numeric(toString(significant[grep(allCorr[i,2], significant[,1]), 2]))  < pval_threshold){
        #remove duplicates
        if(grep(allCorr[i,1], significant[,1]) < grep(allCorr[i,2], significant[,1])){
          #Add Variable1, Variable2, correlation, pvalue of Variable 1, and pvalue of variable 2 to the dataframe
          insert <- cbind(toString(allCorr[i,1]), toString(allCorr[i,2]), as.numeric(toString(allCorr[i,3])), as.numeric(toString(significant[grep(allCorr[i,1], significant[,1]), 2])), as.numeric(toString(significant[grep(allCorr[i,2], significant[,1]), 2])))
          sigCorr <- rbind(sigCorr, insert)
        }
      }
    }
  }
}

#Dendograms##############

#Variable Inputs
tData <- my.data.f.rz
#########################

#Construct Dendograms from various subsets of phenotypes
#find the phenotypes that are in the bettercorr
index=as.numeric(match(as.character(betterCorr[,1]),as.character(names(tData))))
cl.tData=tData[,index]
distance=dist(t(cl.tData))
#Set NA distance to maximum distance
distance[which(is.na(distance))]=max(na.omit(distance))+1
plot(hclust(distance), main="Sex-Correlated Phenotypes, Female, Non-Transformed")

#Easy Correlations#######
cor(cbind(expr.reduced[,1], newpheno$Sex), use = "pairwise.complete.obs")

#ggplot##################
ggplot(my.data, aes(x = TBIL1, y = ct.NEUT2, colour = Diet)) + geom_point()

#naList##################

#Variable Inputs
tData <- oldpheno
#########################

#Find number of NAs in dataset
naList <- data.frame()
for(i in (grep("Coat.Color", colnames(tData))+1):dim(tData)[2]){
  #compare length of the variable to length of variable with NAs omitted
  insert <- cbind(colnames(tData[i]), length(tData[,i]) - length(na.omit(tData[,i])))
  naList <- rbind(naList, insert)
}

#Principle Component Analysis####

#Variable inputs
na.thres <- 100
tData <- oldpheno
pc.thres <- 0.1
#################################

#find all phenotypes in naList with NAs more than na.thres
ls <- subset(naList, as.numeric(as.character(naList[,2])) > na.thres)
#Find indices for each of these phenotypes in data supplied
index=as.numeric(match(as.character(ls[,1]),as.character(names(tData))))
#remove these phenotypes from the dataset
my.data.reduced <- tData[,-index]
#remove any mice that still have NAs from the data
my.data.reduced <- na.omit(my.data.reduced)
#Run PCA using correlation
pc <- princomp(x = my.data.reduced[, c((grep("Coat.Color", colnames(my.data.reduced))+1):dim(my.data.reduced)[2])], cor = TRUE)
#Break up the components and find the base phenotypes
pc.ld <- loadings(pc)
pc.scores <- pc$scores
#separate the first 4 components
comp1 <- pc.ld[,1]
comp2 <- pc.ld[,2]
comp3 <- pc.ld[,3]
comp4 <- pc.ld[,4]
#find all phenotypes in these components with values > 0.1
comp1.pheno <- subset(comp1, abs(comp1) > pc.thres)
comp2.pheno <- subset(comp2, abs(comp2) > pc.thres)
comp3.pheno <- subset(comp3, abs(comp3) > pc.thres)
comp4.pheno <- subset(comp4, abs(comp4) > pc.thres)

#add sex, diet, generation, and litter to scores data
sex.scores <- my.data[as.numeric(rownames(pc.scores)),4]
diet.scores <- my.data[as.numeric(rownames(pc.scores)),7]
gen.scores <- my.data[as.numeric(rownames(pc.scores)),5]
parity.scores <- my.data[as.numeric(rownames(pc.scores)),6]
pc.scores <- cbind(pc.scores, sex.scores, diet.scores, gen.scores, parity.scores)
pc.scores <- pc.scores[,c(68:71, 1:67)]
colnames(pc.scores)[1:4] <- c("Sex", "Diet", "Gen", "parity")

#pc.Significant#################

#Variable Values
pvalue_threshold <- .05
tData <- pc.scores
################################

pc.significant <- data.frame()
for (i in 5:length(colnames(pc.scores))) {
  #Finding linear models for all the variables influence on components
  #1 is Sex, 2 is Diet, 3 is Generation, 4 is Parity
  test1 <- lm(tData[,i]~tData[,1]+tData[,2]+tData[,3] + tData[,4])
  test2 <- lm(tData[,i]~tData[,1]+tData[,3] + tData[,4])
  test3 <- lm(tData[,i]~tData[,2]+tData[,3] + tData[,4])
  test4 <- lm(tData[,i]~tData[,3] + tData[,4])
  test5 <- lm(tData[,i]~tData[,1]*tData[,2] + tData[,1] + tData[,2] + tData[,3] + tData[,4])
  test6 <- lm(tData[,i]~tData[,1]+tData[,2]+tData[,3])
  #tests for na
  if(!is.na(anova(test1, test2)[2,6]) & !is.na(anova(test1, test3)[2,6]) & !is.na(anova(test1, test5)[2,6] < .05)){
    #tests if p-value is significant
    if(anova(test1, test2)[2,6] < pvalue_threshold | anova(test1, test3)[2,6] < pvalue_threshold | anova(test1, test5)[2,6] < pvalue_threshold | anova(test1, test6)[2,6] < pvalue_threshold){
      #create dataframe of p-values for the influence of just sex, just diet, etc. on the phenotype
      insert <- cbind(anova(test1, test2)[2,6], anova(test1, test3)[2,6], anova(test1, test5)[2,6], anova(test1, test6)[2,6], i, colnames(tData)[i])
      pc.significant <- rbind(pc.significant, insert)
    }
  }
}
#name columns
colnames(pc.significant) <- c("Diet Test", "Sex Test", "Inter. Test", "Parity Test", "Index", "Component")

#pcCorr##################

#Variable inputs
phenotype <- "Diet"
corr_threshold <- .4
#########################

pcCorr <- data.frame()
#Find where in the dataframe the testing phenotype is
index <- grep(phenotype, colnames(pc.scores))
for (i in 5:length(colnames(pc.scores))) {
  #find correlations between the phenotype and each component
  cor <- cor(cbind(as.numeric(pc.scores[,index]), pc.scores[,i]), use = "pairwise.complete.obs")
  #Check if NA
  if (!is.na(cor[2,1])){
    #Check if it has a significant p-value and correlation
    if(cor[2,1] > corr_threshold && length(grep(i, pc.significant[,6])) != 0){
      #Create dataframe
      insert <- cbind(colnames(pc.scores)[i], cor[2,1])
      pcCorr <- rbind(pcCorr, insert)
    }
  }
}