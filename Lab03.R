# Lab 03
rm(list=ls())

# Question (1)
data<-c(2,3,2,6,3,5,6,2,6,6,2,6,6,2,3,6,6,6,5,6,6,5,6,6,6,6,6,4,6,3,3,3,6,6,5,6,6)
probPrior <- c(0.01,0.99) 
probLoaded <- c(rep(1/10,5),5/10)
container <- vector()
for (i in 1:length(data)){
  denom <- probPrior[1] * probLoaded[data[i]] + probPrior[2] * 1/6
  probPrior[1] <- probLoaded[data[i]] * probPrior[1]/denom
  container <- c(container,probPrior[1])
  probPrior[2] <- 1 - probPrior[1]
}
plot(1:length(data),container, main= "plot for question #1",
     xlab="number of rolls", ylab="posterior probability that we have picked up a loaded die ", pch=19)

# Question (2)A
probs <- seq(0, 1, 0.001)
maxY <- max(dbeta(probs, 1, 1), dbeta(probs, 6, 6))
minY <- min(dbeta(probs, 1, 1), dbeta(probs, 6, 6))
labX <- "probability"
labY <- "dbeta"
plot(probs, dbeta(probs, 1, 1), ylim=c(minY, maxY), col="red", xlab=labX, ylab=labY)
par(new=TRUE)
plot(probs, dbeta(probs, 6, 6), ylim=c(minY, maxY), col="blue", 
     main= "plot for question #2A", xlab=labX, ylab=labY)
legend("topleft", legend=c("dbeta(probs,1,1)", "dbeta(probs,6,6)"), 
       col=c("red", "blue"), lty=3, lwd = 3, cex=0.8)

# Question (2)B
dbetaUpdateTwoNew_11 <- dbeta(probs,1+1, 1+1)
dbetaUpdateTwoNew_12 <- dbeta(probs,6+1, 6+1)
plot(probs,dbetaUpdateTwoNew_11, xlab=labX, ylab= labY, col="red",
     ylim=c(min(dbetaUpdateTwoNew_11,dbetaUpdateTwoNew_12),
            max(dbetaUpdateTwoNew_11,dbetaUpdateTwoNew_12)))
par(new=TRUE)
plot(probs,dbetaUpdateTwoNew_12, xlab=labX, ylab= labY, col= "blue",
     main= "plot for question #2B: involving the 2 new coin flips",
     ylim=c(min(dbetaUpdateTwoNew_11,dbetaUpdateTwoNew_12), 
            max(dbetaUpdateTwoNew_11,dbetaUpdateTwoNew_12)))
legend("topleft", legend=c("dbeta(probs,1+1,1+1)","dbeta(probs,1+6,1+6)"), 
       col=c("red", "blue"), lty=2, lwd = 3, cex=0.5)

dbetaUpdateTwoNew_21 <- dbeta(probs,1+400, 1+400)
dbetaUpdateTwoNew_22 <- dbeta(probs,6+400, 6+400)
plot(probs,dbetaUpdateTwoNew_21, col="red", pch=21, xlab= labX, ylab= labY,
     ylim=c(min(dbetaUpdateTwoNew_21,dbetaUpdateTwoNew_22),
            max(dbetaUpdateTwoNew_21,dbetaUpdateTwoNew_22)))
par(new=TRUE)
plot(probs,dbetaUpdateTwoNew_22, col="blue", pch=20, xlab= labX, ylab =labY,
     main= "plot for question #2B: involving the 800 new coin flips",
     ylim=c(min(dbetaUpdateTwoNew_21,dbetaUpdateTwoNew_22),
            max(dbetaUpdateTwoNew_21,dbetaUpdateTwoNew_22)))
legend("topleft", legend=c("dbeta(probs,1+400,1+400)","dbeta(probs,6+400, 6+400)"), 
       col=c("red", "blue"), lty=2, lwd = 3, cex=0.5)
##########################################################################################
