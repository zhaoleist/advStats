# Lab 4
rm(list=ls())

# Question (1A)
probs <- seq(0,1,0.01)
aNum <- 0.9932621
expDist <- dexp(probs, rate=5)/aNum
plot(probs, expDist, ylab="density", main="question 1(A)")

####################################################################################
# define metropolis function
metropolis_priorExp <- function(piOld=0.5, iteration=100000, numHeads, numTrials){
  container <- vector()
  aNum <- 0.9932621
  for (i in 1:iteration){
    pOld <- dexp(piOld, rate=5)/aNum * dbinom(numHeads, numTrials, piOld)
    piNew <- piOld + rnorm(1, 0, sd=0.01)
    if (piNew > 1)
      piNew <- 1
    if (piNew < 0)
      piNew <- 0
    pNew <- dexp(piNew, rate=5)/aNum * dbinom(numHeads, numTrials, piNew)
    ratio <- pNew/pOld
    if (ratio > 1 || ratio >= runif(1))
      piOld <- piNew
    container[i] <- piOld
  }
  return(container)
}

#####################################################################################
# define grid function
grid_priorExp <- function(xVals, i=1, mySum=0, numHeads, numTrials){
  
  posteriorDist <- vector()

  for (x in xVals){
    posteriorDist[i] <- dexp(x, rate=5)/aNum * dbinom(numHeads, numTrials, x)
    mySum <- mySum + posteriorDist[i]
    i <- i + 1
  }
  return(c(posteriorDist, mySum))
}

###################################################################################

# Question (1B)
posteriorDist_metro <- metropolis_priorExp(0.5, 100000,14,24)
myHist <- hist(posteriorDist_metro, breaks=200, plot=FALSE)

myVector <- grid_priorExp(myHist$mids,1,0,14,24)
mySum <- myVector[length(myVector)]
posteriorDist_grid <- myVector[-length(myVector)]

minY <- 0
maxY <- max(myHist$counts/sum(myHist$counts), posteriorDist_grid/mySum, dbeta(myHist$mids, 40+14, 40+10)/sum(dbeta(myHist$mids, 40+14, 40+10)))

plot(myHist$mids, myHist$counts/sum(myHist$counts), xlim=c(0,1), col="green", ylim=c(minY, maxY),
     main="with additional observation of 14 heads & 10 tails ")
points(myHist$mids, posteriorDist_grid/mySum, xlim=c(0, 1),  col="blue", xlab="", ylab="", ylim=c(minY, maxY))
points(myHist$mids, dbeta(myHist$mids, 40+14, 40+10)/sum(dbeta(myHist$mids, 40+14, 40+10)), 
     xlim=c(0,1), ylim=c(minY, maxY), col="red", xlab="", ylab="")
legend("topleft", legend=c("metropolis","grid", "prior dbeta(40, 40)"), 
       col=c("green", "blue", "red"), lty=3, lwd = 3, cex=0.6)
######################################################################################

# Question (1C)
posteriorDist_metro_2 <- metropolis_priorExp(0.5, 500000,583,1000)
myHist_2 <- hist(posteriorDist_metro_2, breaks=200, plot=FALSE)
myVector_2 <- grid_priorExp(myHist_2$mids,1,0,583,1000)
mySum_2 <- myVector_2[length(myVector_2)]
posteriorDist_grid_2 <- myVector_2[-length(myVector_2)]
minY_2 <- 0
maxY_2 <- max(myHist_2$counts/sum(myHist_2$counts), posteriorDist_grid_2/mySum_2, dbeta(myHist_2$mids, 40+583, 40+417)/sum(dbeta(myHist_2$mids, 40+583, 40+417)))

plot(myHist_2$mids, myHist_2$counts/sum(myHist_2$counts), xlim=c(0.4,0.8), ylim=c(minY_2, maxY_2), col="green",
     main="with additional observation of 583 heads & 417 tails")

points(myHist_2$mids,posteriorDist_grid_2/mySum_2, xlab="", ylab="", col="blue", xlim=c(0.4,0.8), ylim=c(minY_2, maxY_2))

points(myHist_2$mids, dbeta(myHist_2$mids, 40+583, 40+417)/sum(dbeta(myHist_2$mids, 40+583, 40+417)), 
     xlim=c(0.4,0.8), ylim=c(minY_2, maxY_2), col="red", xlab="", ylab="")
legend("topleft", legend=c("metropolis", "grid", "prior dbeta(40, 40)"), 
       col=c("green", "blue", "red"), lty=3, lwd = 3, cex=0.6)