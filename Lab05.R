# lab 5
rm(list=ls())

# (1)
myT <- read.table("nc101_scaff_dataCounts.txt", header=TRUE, row.names = 1)

# (2)
par(mfrow=c(3,2))
myT_log10 <- log10(myT+1)
plot(myT_log10[,1], myT_log10[,2], main="log scale counts for the two samples")
lines(c(0,4), c(0,4), col="red")
# ANSWER:
# The biological replicates do not seem to have similar patterns of gene expression. 
# While the red line indicates the case that the two replicates have same patterns, most actual points are to the right of the line, which means most genes 
# have more expression in replicate 1.

# (3) (4)
sumCol <- apply(myT, 2, sum)
pValues <- vector(length=nrow(myT))
for (i in 1:nrow(myT)){
  m <- matrix(c(myT[i,1], myT[i, 2], sumCol[1]-myT[i,1], sumCol[2]-myT[i,2]), nrow=2, byrow=TRUE)
  pValues[i] <- fisher.test(m)$p.value
}
hist(pValues, breaks=20, main="p-values for fisher.test")
# The p-value for question #(3) is 1.670017e-11.

myT_rmvLowAbund <- myT[myT$D2_01 + myT$D2_02 > 50,]
sumCol_2 <- apply(myT_rmvLowAbund, 2, sum)
pValues_2 <- vector(length=nrow(myT_rmvLowAbund))
for (i in 1:nrow(myT_rmvLowAbund)){
  m <- matrix(c(myT_rmvLowAbund[i,1], myT_rmvLowAbund[i, 2], sumCol_2[1]-myT[i,1], sumCol_2[2]-myT[i,2]), nrow=2, byrow=TRUE)
  pValues_2[i] <- fisher.test(m)$p.value
}
hist(pValues_2, breaks=20, main="p-values for fisher.test without low abundance counts")
# Answer for #(4):
# They are not uniformly distributed and I don't expect them to be. If we say p <- 0.05 is "significant", 
# the ratio of significant p-values is: 29.1% ( hist(pValues, breaks=20, main="p-values for fisher.test")$counts[1]/length(pValues) ), 
# so we say they are less significant. After removing low abundance expression, 1) the histogram lost the right-most bar, which was caused 
# by the raw zero count in the dataset; 2) the ratio of significant p-values went up to 
# 43.9% (hist(pValues_2, breaks=20, main="p-values for fisher.test")$counts[1]/length(pValues_2)). 

# (5) (6)
myT <- myT + 1
freq_expt <- myT[1,1]/sum(myT[,1])
poisson.test(myT[1,2], sum(myT[,2]), freq_expt)$p.value
# The p-value for question #(5) is: 1.139341e-13
sumCol_plusOne <- apply(myT, 2, sum)
pValues_3 <- vector(length=nrow(myT))
for (i in 1:nrow(myT)){
  freq_background <- myT[i,1]/sumCol_plusOne[1]
  pValues_3[i] <- poisson.test(myT[i,2], sumCol_plusOne[2], freq_background)$p.value
}
hist(pValues_3, breaks=20, main="p-values for possion.test")
plot(pValues, pValues_3, main="p-values from fish.test vs. p-values from possion.test")
lines(c(0,1), c(0,1), col="red")
# Answer for #(6):
# The histograms look similar, but from the last graph, they don't agree.
