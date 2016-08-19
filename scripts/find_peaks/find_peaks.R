#Load in qtl scan data

load("scans.RData")
load("Churchill_project/Results/diffs.Rdata")
load("allcoef.Rdata")
load("med_scans.Rdata")

#save qtl scans for each phenotype to a pdf
#5 scans per graphics window



pdf(paste("Churchill_project/results/interaction_graphs/", gene, "_interaction_scans.pdf"))
for(i in 1:128){
  gene <- colnames(allscan1.add$lod)[i]
  par(mfrow = c(3, 2), mar = c(1,2,2,1))
  plot(allscan1.intsex[,i], main = paste(gene, " with intcovar sex"))
  plot(allscan1.intdiet[,i], main = paste(gene, " with intcovar diet"))
  plot(allscan1.difsex[,i], main = paste(gene, " difference"))
  plot(allscan1.difdiet[,i], main = paste(gene, " difference"))
  plot(allscan1.add[,i], main = paste(gene, " additive"))
  #graphics.off()
}

#find peaks for the qtl scans
#additive model uses a threshold of 7 because of previous knowledge
#perms on interactive model showed that the 95% significance threshold was around 10, 
  #but we lowered the threshold to 9 to get anything that was close
#perms on difference between the models showed the 95% threshold to be around 6, but
  #we lowered the threshold to 5 for the same reason as above
add <- find_peaks(allscan1.add, threshold = 7)
intsex <- find_peaks(allscan1.intsex, threshold = 9)
intdiet <- find_peaks(allscan1.intdiet, threshold = 9)
difsex <- find_peaks(allscan1.difsex, threshold = 5)
difdiet <- find_peaks(allscan1.difdiet, threshold = 5)

#finds all significant peaks for intsex that have a corresponding significant difference peak
intpeaks.sex <- data.frame()
x <- 1
while(x < (nrow(intsex) + 1)){
  #creates a list of each instance of the phenotype in intsex
  index <- grep(intsex[x,2], intsex[,2])
  #compares each instance of the phenotype in difsex to each in intsex
  for(i in grep(intsex[x,2], difsex[,2])){
    for(j in index){
      #checks if the peak is on the same chromosome and near the same position
      if((as.character(intsex[j,3]) == as.character(difsex[i,3])) & as.numeric(intsex[j,4]) < 
           (as.numeric(difsex[i,4]) + 5) & 
           as.numeric(intsex[j,4]) > (as.numeric(difsex[i,4]) - 5)){
        insert <- cbind(intsex[j,2], intsex[j,3])
        intpeaks.sex <- rbind(intpeaks.sex, insert)
      }
    }
  }
  #moves to next phenotype
  x <- index[length(index)] + 1
}

#Same as above for intdiet
intpeaks.diet <- data.frame()
x <- 1
while(x < (nrow(intdiet) + 1)){
  index <- grep(intdiet[x,2], intdiet[,2])
  for(i in grep(intdiet[x,2], difdiet[,2])){
    for(j in index){
      if((as.character(intdiet[j,3]) == as.character(difdiet[i,3])) & as.numeric(intdiet[j,4]) < 
         (as.numeric(difdiet[i,4]) + 5) & 
         as.numeric(intdiet[j,4]) > (as.numeric(difdiet[i,4]) - 5)){
        insert <- cbind(intdiet[j,2], intdiet[j,3])
        intpeaks.diet <- rbind(intpeaks.diet, insert)
      }
    }
  }
  x <- index[length(index)] + 1
}

#find peaks that are in both intsex and intdiet 
intpeaks.both <- data.frame()
for(i in 1:nrow(intpeaks.sex)){
  index <- grep(intpeaks.sex[i,1], intpeaks.diet[,1])
  #commented code allows to only check for peaks at the same position
  if(length(index) != 0){# & length(grep(intpeaks.sex[i,2], intpeaks.diet[index,2])) != 0){
    insert <- cbind(as.character(intpeaks.sex[i,1]))
    intpeaks.both <- rbind(intpeaks.both, insert)
  }
}

intpeaks.noadd.samepos <- data.frame()
for(i in 1:nrow(intsex)){
  index.diet <- grep(intsex[i,2], intdiet[,2])
  index.add <- grep(intsex[i,2], add[,2])
  for(j in index.diet){
    if(intsex[i,3] == intdiet[j,3] & 
       as.numeric(intdiet[j,4]) < (as.numeric(intsex[i,4]) + 5) & 
       as.numeric(intdiet[j,4]) > (as.numeric(intsex[i,4]) - 5) &
       length(grep(intsex[i,3], add[index.add,3])) == 0){
      insert <- cbind(intsex[i,2], intsex[i,3])
      intpeaks.noadd.samepos <- rbind(intpeaks.noadd.samepos, insert)
    }
  }
}

intpeaks.noadd <- data.frame()
for(i in 1:nrow(intsex)){
  index.diet <- grep(intsex[i,2], intdiet[,2])
  index.add <- grep(intsex[i,2], add[,2])
  for(j in index.diet){
    if(intsex[i,3] == intdiet[j,3] & 
       length(grep(intsex[i,3], add[index.add,3])) == 0){
      insert <- cbind(intsex[i,2], intsex[i,3])
      intpeaks.noadd <- rbind(intpeaks.noadd, insert)
    }
  }
}

coef_jumps_intsex <- data.frame()
for(i in 1:length(allcoef.intsex)){
  for (j in 1:19) {
    point <- bayes_int(allscan1.intsex[,i], chr = j)[1,2]
    marker <- find_marker(allscan1.intsex$map, chr = j, pos = point)
    sort <- allcoef.intsex[[i]][[j]]$coef[marker,13:19]
    sort <- sort(sort)
    for(k in 1:6){
      dif <- sort[k+1] - sort[k]
      if(dif > 1){
        insert <- cbind(names(allcoef.add)[i], j, dif)
        coef_jumps_intsex <- rbind(coef_jumps_intsex, insert)
      }
    }
  }
}

coef_jumps_intdiet <- data.frame()
for(i in 1:length(allcoef.intdiet)){
  for (j in 1:19) {
    point <- bayes_int(allscan1.intdiet[,i], chr = j)[1,2]
    marker <- find_marker(allscan1.intdiet$map, chr = j, pos = point)
    sort <- allcoef.intdiet[[i]][[j]]$coef[marker,13:19]
    sort <- sort(sort)
    for(k in 1:6){
      dif <- sort[k+1] - sort[k]
      if(dif > 1){
        insert <- cbind(names(allcoef.add)[i], j, dif)
        coef_jumps_intdiet <- rbind(coef_jumps_intdiet, insert)
      }
    }
  }
}

drops.intsex <- data.frame()
for (i in 2) {
  for (j in 1:128) {
    for (k in 1:19) {
      med_chr <- med_scans[[i]][[j]][[k]][which(med_scans[[i]][[j]][[k]][,3] == as.character(k)),]
      dif <- mean(med_scans[[i]][[j]][[k]][,9]) - min(med_chr[,9])
      if(dif > 3.077478){
        insert<- cbind(names(med_scans)[i], names(med_scans[[i]])[j], 
             names(med_scans[[i]][[j]])[k], dif, 
             med_chr[grep(min(med_chr[,9]), med_chr[,9]),2])
        drops.intsex <- rbind(drops.intsex, insert)
      }
    }
  }
}
