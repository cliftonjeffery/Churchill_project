library(qtl2scan)
library(qtl2plot)
library(qtl2geno)
setwd("~/Desktop/Rstudio")
#read_cross2 to load in data
write_control_file("DO_control_file.json", crosstype = "do", 
                   geno_file = "expr.csv", #"Churchill_project/Data/annot.csv", 
                   pheno_file = "Churchill_project/Data/DO_phenotype_data.csv", 
                   gmap_file = "Churchill_project/Data/DO_gmap.csv",
                   covar_file = "Churchill_project/Data/DO_covar.csv", 
                   sex_file = "Churchill_project/Data/DO_covar_sex.csv", 
                   sex_codes = c("1"="female", "0"="male"))
load("Churchill_project/Data/DO478_ExprData_4eQTLMapping.Rdata")
covar.data <- read.csv("Churchill_project/Data/DO_covar.csv", header = T)
phenocovar.data <- read.csv("Churchill_project/Data/DO_phenocovar.csv", header = T)
geno.data <- read.csv("Churchill_project/Data/geno_file.csv", header = T)
covar.data$Sex <- gsub("F", toString(1), covar.data$Sex) 
covar.data$Sex <- gsub("M", toString(0), covar.data$Sex)
phenocovar.data$Sex <- gsub("F", toString(1), phenocovar.data$Sex) 
phenocovar.data$Sex <- gsub("M", toString(0), phenocovar.data$Sex)

DO <- read_cross2("DO_control_file.json")
#calculate the probabilitiy of conditional genotypes
#does not insert back into the cross object, but returned as a 3D array
pr <- calc_genoprob(DO, step=1, err=0.002)
#calculate kinship of mice (allows choice of kinship methods)
kinship <- calc_kinship(probs.478)
#get covariates for null hypothesis when performing QTL on the X chromosome
Xcovar <- get_x_covar(iron)
#Genome scan using Haley-Knott (kinship can also be included, as well as other covars)
scan <- scan1(pr, iron$pheno, Xcovar = Xcovar)
quartz()
#parameters of border (1 is bottom, 2 is left side, 3 is top, and 4 is right side)
par(mar=c(5.1, 4.1, 1.1, 1.1))
#find maximum LOD score
ymx <- maxlod(scan) # overall maximum LOD score
#plot the first column of the scan (only one column can be graphed at a time)
plot(scan, lodcolumn=1, col="slateblue", ylim=c(0, ymx*1.02))
#plot the second column, adding on column 2
plot(scan, lodcolumn=2, col="violetred", add=TRUE)
#set location and contents of legend
legend("topleft", lwd=2, col=c("slateblue", "violetred"), colnames(scan$lod), bg="grey90")
#find peaks on each chromosome
#can find two peaks on one chromosome with peakdrop =, but apparently you're supposed to be carful with that
find_peaks(scan, threshold = 4, peakdrop = 1.8, drop = 1.5)
#find area where the QTL is 95% sure to be located
#again can only find two peaks on same chromosome with peakdrop
bayes_int(scan, lodcolumn=1, chr=7, peakdrop = 1.8, prob=0.95)

#Shows QTL effects
c2eff <- scan1coef(pr[,"2"], iron$pheno[,"liver"])
col <- c("slateblue", "violetred", "green3")
quartz()
par(mar=c(4.1, 4.1, 1.1, 2.6), las=1)
plot(c2eff, columns = 1:3)
co <- c2eff$coef
for(i in 1:ncol(co))
  axis(side=4, at=co[nrow(co),i], colnames(co)[i], tick=FALSE, col.axis=col[i])

#load new data (DO)
file <- paste0("https://raw.githubusercontent.com/rqtl/", "qtl2data/master/DOex/DOex.zip")
DOex <- read_cross2(file)
#again run genoprobs
pr <- calc_genoprob(DOex, error_prob=0.002)
#find allele probs
apr <- genoprob_to_alleleprob(pr)
#again find kinship
k <- calc_kinship(apr, "loco")
#Turn sex into a covariate
sex <- (DOex$covar$Sex == "male")*1
names(sex) <- rownames(DOex$covar)
#perform the QTL scan
out <- scan1(apr, DOex$pheno, k, sex)

#plot lod peaks
par(mar=c(4.1, 4.1, 0.6, 0.6))
plot(out)

#plot contributions of each parent strain
coef_c2 <- scan1coef(apr[,"2"], DOex$pheno, k[["2"]], sex)
quartz()
par(mar=c(4.1, 4.1, 0.6, 0.6))
plot_coefCC(coef_c2, bgcolor="gray95")
legend("bottomright", col=CCcolors, names(CCcolors), ncol=2, lwd=2, bg="gray95")

#Find the locatiion of the marker under the peak
marker <- rownames(max(out, chr="2"))
peak_Mbp <- DOex$pmap[["2"]][marker]

#download SNP data
tmpfile <- tempfile()
file <- "https://raw.githubusercontent.com/rqtl/qtl2data/master/DOex/c2_snpinfo.rds"
download.file(file, tmpfile, quiet=TRUE)
snpinfo <- readRDS(tmpfile)
unlink(tmpfile)

#plot lod score of each QTL
snpinfo$sdp <- calc_sdp(snpinfo[,-(1:4)])
apr_Mbp <- apr
apr_Mbp$map <- DOex$pmap
snppr <- genoprob_to_snpprob(apr_Mbp, snpinfo)
out_snps <- scan1(snppr, DOex$pheno, k[["2"]], sex)
quartz()
par(mar=c(4.1, 4.1, 0.6, 0.6))
plot_snpasso(out_snps)
#output all QTLs with LOD scores less than drop lower than the maximum score
top_snps(out_snps, drop = .5)
#highlight all the QTLs on the graph with LOD scores less than frop lower than the maximum score
plot_snpasso(out_snps, drop = 1.5)