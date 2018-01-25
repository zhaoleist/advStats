rm(list=ls())

# Question 1 with formulation # 2
sampleSize <- 30
probMut <- 1/3
dbinom(12, sampleSize, probMut)    # 0.01101246
plot(0:sampleSize, dbinom(0:sampleSize, sampleSize, probMut), xlab=" observed number of people with the mutation ", ylab= "density")
myMean <- sampleSize * probMut    
myVars <- sampleSize * probMut * (1 - probMut)    
#ANSWER:
" The probility of seeing exatly 12 people with the mutation is 0.01101246.
  Mean for the expected number of people with the mutation is 10.
  Variance for the expected number of people with the mutation is 6.666667. "

# 3A
set.seed(1000)
numPatients <- 10000
numExperiments <- 1000
myVals <- rbinom(numExperiments,numPatients,0.5)    # one liner: rbinom(1000,10000,0.5)

# 3B
expectedMean <- numPatients * 0.5
expectedVars <- numPatients * 0.5 * (1 - 0.5)
actualMean <- mean(myVals)
actualVars <- var(myVals)
# ANSWER
" The expected mean is 5000 and the expected variance is 2500.
  The actuial mean is 5000.688 and the actual variance is 2700.287."

# 3C
pValues <- vector(length=length(myVals), mode="double")
for (i in 1:length(myVals)){
    pValues[i] <- binom.test(myVals[i], numPatients)$p.value
}
hist(pValues)
# ANSWER
"I expected to see uniform distribution and 
the p-values of the histogram are approximate uniform distribution."

# 3D
newProb <- 2/3
pValues_2 <- vector(length=length(myVals), mode="double")
for (p in 1:length(myVals)){
  pValues_2[p] <- binom.test(myVals[p], numPatients, p= newProb)$p.value
}
hist(pValues_2)    # with p = 2/3

pValues_3 <- vector(length=length(myVals), mode="double")
for (q in 1:length(myVals)){
  pValues_3[q] <- binom.test(myVals[q], numPatients, p=0.49)$p.value
}
histo(pValues_3)    # with p = 0.49
# ANSWER
" Changing the expected value from 1/2 to 2/3, all the p-values become very small. Interestingly,
  if we change 1/2 to a much bigger value, like 3/4, we end up with the same extremely small 
  p-values in R.
  
  If the expected value changed to 0.49 or 0.51, I don't expect to see the same shape of the
  p-value histogram because the actual distribution is no longer symmatric around the distribution 
  with the expected value of 0.49 or 0.51. By histogram, we can verify that we see more of 
  small p-values."

