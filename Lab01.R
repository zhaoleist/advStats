# a "loaded" die that has a 10% chance of getting a 1-5 and a 50% chance of getting a 6.
rm(list=ls())

# 1 What is the mean and variance for the loaded dice?
loaded <- 1:6
myMean <- sum(1/10 * (1:5)) + 1/2 * 6
myVar <- sum((loaded[1:5] - myMean)^2) * 1/10 + 1/2 *(6 - myMean)^2
# ANSWER:
" mean: 4.5, variance: 3.25"

# 2 Make a function in R that "rolls" this dice; return a vector containing the rolls. 
myRolls <- function(times){
  outcomes <- vector(length=times, mode="double")
  for (i in 1: times){
    outcomes[i] <- sample(1:6, 1, prob = c(rep(1/10,5),1/2))
  }
  return(outcomes)
}
iter <- 10000
myRolls(iter)

# 3 Make a histogram of some large number of rolls.  
# Do the rolls of the loaded die approximate a uniform distribution?
hist(myRolls(iter), breaks=0:6)
barplot(table(myRolls(iter)), xlab="roll", ylab="frequency")
# Since the rolls are discrete, it's more appropriate to use barplot()
# ANSWER:
" The rolls of the loaded die does not approximate a uniform distribution. We see more of a 6"

# 4 number of rolls for convergence
trialSizes <- c(5,10,15,20,25,30,40,50,100,200,300,400,500,1000,2000,3000,4000,5000,
                10000,20000,30000,40000,50000,100000)
means <- c()
vars <- c()
for (m in trialSizes){
  rolls <- vector(length=m, mode="double")
  rolls <- myRolls(m)
  means <- c(means,mean(rolls))
  vars <- c(vars, var(rolls))
}
myLabels <- factor(1:length(trialSizes))
plot(log10(trialSizes), means)
text(log10(trialSizes), means, label=myLabels, pos=1, cex=0.5, col="red")
lines(log10(trialSizes),rep(myMean,length(trialSizes)))
windows()
plot(log10(trialSizes), vars)
text(log10(trialSizes), vars, label=myLabels, pos=1, cex=0.5, col="red")
lines(log10(trialSizes),rep(myVar,length(trialSizes)))
# ANSWER:
" After a few attemps, it seems that both means vs. trial size plot and 
 variance vs. trial size plot started getting convergence on the expected values
 from the 15th value of the trialSize vector, which is 2000." 

####################################################################################