library(qtl2geno)
library(qtl2plot)
library(qtl2scan)
library(ggplot2)
library(intermediate)

setwd("Desktop/Rstudio")
load("allcoef.Rdata")
load("Churchill_project/results/scans.Rdata")
load("Churchill_project/results/diffs.Rdata")
load("Churchill_project/results/med_scans.Rdata")

#create a pdf with 4 graphs per window of coefficient scans and their corresponding qtl graphs
for(i in 1:length(allcoef.add)){
  #creat a new pdf for each phenotype
  pdf(paste("Churchill_project/Graphs/coefscans/", names(allcoef.add)[i], ".pdf", sep = ""))
  #graph each chromosome individually
  for(j in 1:19){
    par(mfrow = c(3,2), mar = c(2,2,1,1))
    plot_coef(allcoef.intsex[[i]][[j]][,13:19], col = CCcolors[2:8], main = 
                  paste("intcoefs for sex interactive model for ", names(allcoef.add)[i], sep = ""))
    point <- bayes_int(allscan1.intsex[,i], chr = j)[1,2]
    abline(v = point, lwd = 2, lty = 2)
    plot_coef(allcoef.intdiet[[i]][[j]][,13:19], col = CCcolors[2:8], main = 
                paste("intcoefs for diet interactive model for ", names(allcoef.add)[i], sep = ""))
    point <- bayes_int(allscan1.intdiet[,i], chr = j)[1,2]
    abline(v = point, lwd = 2, lty = 2)
    plot_coef(allcoef.intsex[[i]][[j]][,1:8], col = CCcolors, main = 
                  paste("addcoefs for sex interactive model for ", names(allcoef.add)[i], sep = ""))
    #Plot a line on the coefficient scan at the qtl peak for that model
    point <- bayes_int(allscan1.intsex[,i], chr = j)[1,2]
    abline(v = point, lwd = 2, lty = 2)
    plot_coef(allcoef.intdiet[[i]][[j]][,1:8], col = CCcolors, main = 
                  paste("addcoefs for diet interactive model for ", names(allcoef.add)[i], sep = ""))
    point <- bayes_int(allscan1.intdiet[,i], chr = j)[1,2]
    abline(v = point, lwd = 2, lty = 2)
    plot_coefCC(allcoef.add[[i]][[j]], main = 
                  paste("additive model for ", names(allcoef.add)[i], sep = ""))
    point <- bayes_int(allscan1.add[,i], chr = j)[1,2]
    abline(v = point, lwd = 2, lty = 2)
    thr <- max(c(find_peaks(allscan1.add[,i], threshold = 0)[j,5],
                 find_peaks(allscan1.intsex[,i], threshold = 0)[j,5], 
                 find_peaks(allscan1.intdiet[,i], threshold = 0)[j,5]), na.rm = TRUE)
    #plot all three qtl graphs in the same window
    plot(allscan1.add[,i], ylim = c(0, (thr + 1)), chr = j, col = "violetred", 
                 main = paste("chr", j))
    plot(allscan1.intdiet[,i], ylim = c(0, (thr + 1)), chr = j, col = "deepskyblue", add = TRUE)
    plot(allscan1.intsex[,i], ylim = c(0, (thr + 1)), chr = j, col = "forestgreen", add = TRUE)
    legend("topright", col = c("violetred", "deepskyblue", "forestgreen"), c("additive", "intdiet", "intsex"), lwd = 2)
  }
  graphics.off()
}

intpeaks <- intpeaks.diet
pdf("diet_interactive_mediation_scans.pdf")
for(i in 1:nrow(intpeaks)){
  par(mfrow = c(3,1))
  plot(allscan1.add[,intpeaks[i,1]], main = paste("additive model for", intpeaks[i,1]))
  plot(allscan1.intdiet[,as.character(intpeaks[i,1])], main = paste("diet interactive model for", intpeaks[i,1]))
  chr <- as.numeric(as.character(intpeaks[i,2]))
  index <- grep(intpeaks[i,1], names(med_scans[[2]]))
  plot(med_scans[[2]][[index]][[chr]])
}
graphics.off()















