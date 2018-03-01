# lab 6
rm(list=ls())
library(DESeq)
myT <- read.table("nc101_scaff_dataCounts.txt",header=TRUE,row.names=1)
# file has inconsistent upper/lowercase column names for W/w. may cause trouble later. Change them to lowercase
names(myT)[4:11] <- tolower(names(myT)[4:11])    # Uppercase W to lowercase w    
par(mfrow=c(2,2))
# normalization by DESeq
timepoints <- c("D2","D2","D2", "w12","w12","w12","w20","w20","w20","w20","w20")
cds <- newCountDataSet(myT, timepoints)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds,sharingMode="gene-est-only")
myTNorm <- counts(cds,normalized=TRUE) + 1

# (5)
vars <- apply(myTNorm, 1, var)
means <- apply(myTNorm, 1, mean)
plot(log10(vars), log10(means), main="log scale variance vs. means")
lines(c(0,6), c(0,6), col="red")
# ANSWER:
# Means do not equal variance and means are smaller than variance.

# (6)
col_2_weeks <- 1:3
col_20_weeks <- 7:11
numAverageSeqs_2_weeks <- sum(apply(myTNorm[,col_2_weeks], 2, sum))/length(col_2_weeks)
numAverageSeqs_20_weeks <- sum(apply(myTNorm[,col_20_weeks], 2, sum))/length(col_20_weeks)

pValuesPoisson <- vector()
averageNum_20_weeks <- vector()
averageNum_2_weeks <- vector()
for (i in 1: nrow(myTNorm)){
  frac_2_weeks <- mean(myTNorm[i,col_2_weeks])/numAverageSeqs_2_weeks
  mean_20_weeks <- mean(myTNorm[i, col_20_weeks])
  averageNum_20_weeks[i] <- mean_20_weeks
  averageNum_2_weeks[i] <- mean(myTNorm[i, col_2_weeks])
  pValuesPoisson[i] <- poisson.test(round(mean_20_weeks), numAverageSeqs_20_weeks, frac_2_weeks, 
                          alternative="two.sided")$p.value

}
plot(log10(averageNum_2_weeks), log10(averageNum_20_weeks), 
     col=ifelse(p.adjust(pValuesPoisson,method="BH") < 0.1, "red","black" ), 
     xlab="log10(average # of reads at two weeks)", ylab="log10(avarage # of reads at twenty weeks)",
     main="from poisson")
lines(c(0,10),c(0,10), col="blue")

# (7)
cds <- newCountDataSet(myT, timepoints)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds,sharingMode="gene-est-only")
res <- nbinomTest(cds, "D2", "w20")
means <- apply(counts(cds,normalized=TRUE), 1,mean)
perGeneEstimates <- fitInfo(cds)$perGeneDispEsts * means * means

# check indices for NA and where mean greater than variance
if (all(!is.nan(perGeneEstimates))==FALSE){
  indices <- which(is.nan(perGeneEstimates))
  perGeneEstimates[indices] <- averageNum_2_weeks[indices] + 1e-08} 
  
indices_2 <- which(perGeneEstimates < averageNum_2_weeks)
perGeneEstimates[indices_2] <- averageNum_2_weeks[indices_2] + 1e-08

pValuesNegativeBinomial <- vector()
myP <- vector()
myR <- vector()
for (i in 1: length(perGeneEstimates)){
  myP[i] <- averageNum_2_weeks[i]/perGeneEstimates[i]
  myR[i] <- averageNum_2_weeks[i]^2/(perGeneEstimates[i]-averageNum_2_weeks[i])
  pValue <- pnbinom(round(averageNum_20_weeks[i]), myR[i], myP[i])
  if (pValue <= 0.5)
    pValuesNegativeBinomial[i] <- pValue * 2
  else if (pValue > 0.5)
    pValuesNegativeBinomial[i] <- (1- pValue) * 2 
}
plot(log10(averageNum_2_weeks), log10(averageNum_20_weeks), 
     col=ifelse(p.adjust(pValuesNegativeBinomial,method="BH") < 0.1, "red","black"),
     xlab="log10(average # of reads at two weeks)", ylab="log10(avarage # of reads at twenty weeks)",
     main="from negative binomial")
lines(c(0,10),c(0,10), col="blue")
# ANSWER: 
# Negative binomial test is more conservative and more appropriate. The variance are greater than means in rna-seq data, 
# which violated the assumption of poisson distribution that mean equals variance. If we use poisson test, variance are
# under-represented and we would end up with more significant p-values.

# (8)
plot(-log10(res$pval), -log10(pValuesNegativeBinomial))
lines(c(0,50), c(0,50), col="blue")
