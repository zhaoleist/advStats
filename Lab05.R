# lab 5
rm(list=ls())

# (1)
myT <- read.table("nc101_scaff_dataCounts.txt", header=TRUE, row.names = 1)
par(mfrow=c(2,2))
# (2)
myT_log10 <- log10(myT+1)
plot(myT_log10[,1], myT_log10[,2], col=c("black","red"), main="Counts for the two samples")
legend("topleft", legend=c(names(myT)[1], names(myT)[2]), 
       col=c("black", "red"), cex=1, pch=21)
# ANSWER:
# Qualitativly, majority of the genes seem to have similar pattern. But we can see more "zero" expression 
# in D2_02 and some more "high" expression in D2_01. 


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
# They are not uniformly distributed. I don't expect them to be. If we say p <- 0.05 as "significant", 
# the ration of significant p-values is: 43.9%, so we say they are less significant. After remove low 
# abundance expression, the histogram lost the right-most bar, which was cause by zero in the dataset.

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
# Answer for #(6):
# They agree.
