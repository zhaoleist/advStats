# lab 11
rm(list=ls())

# (1)
path <- "http://afodor.github.io/classes/stats2015/prePostPhylum.txt"
myT <- read.table(path, header=TRUE, sep="\t")
myTData <- myT[,5:10] 
myPCOA <- princomp(myTData)

# (2)
varExplained <- myPCOA$sdev^2/sum(myPCOA$sdev^2)
pdf("Lab11_Plots_Lei.pdf")

plot(myPCOA$scores[,1], myPCOA$scores[,2], col=factor(myT$genotype), pch=19, main="Colored by genotype", 
     xlab=paste0("PC 1 (", format(varExplained[1] * 100, digits=4), "%)"), 
     ylab=paste0("PC 2 (", format(varExplained[2] * 100, digits=4), "%)"))
legend(-1.3,2, pch=19, legend=levels(factor(myT$genotype)), col=1:2,cex =1, bg="lightblue" )

plot(myPCOA$scores[,1], myPCOA$scores[,2], col=factor(myT$cage), pch=19, main="Colored by cage", 
     xlab=paste0("PC 1 (", format(varExplained[1] * 100,digits=4), "%)"), 
     ylab=paste0("PC 2 (", format(varExplained[2] * 100,digits=4), "%)"))
legend(-1.7,2, pch=19, legend=levels(factor(myT$cage)), col=1:length(unique(myT$cage)),cex = 0.8)

plot(myPCOA$scores[,1], myPCOA$scores[,2], col=factor(myT$time), pch=19, main="Colored by time", 
     xlab=paste0("PC 1 (", format(varExplained[1] * 100,digits=4), "%)"), 
     ylab=paste0("PC 2 (", format(varExplained[2] * 100,digits=4), "%)"))
legend(-1.3,2, pch=19, legend=levels(factor(myT$time)), col=1:length(unique(myT$time)),cex = 1, bg="lightblue")

dev.off()

# (3)
df <- data.frame(matrix(nrow=3, ncol=2))
names(df) <- c("PCA1","PCA2")
row.names(df) <- c("Cage", "Genotype", "Time(pre vs. post)")
indexWT <- which(myT$genotype == "WT" )
index10 <- which(myT$genotype == "10-/-")
indexPre <- which(myT$time == "PRE" )
indexPost <- which(myT$time == "POST")
for (i in 1:ncol(df)){
  myLm <- anova(lm(myPCOA$scores[,i] ~ myT$cage))
  df[1,i] <- myLm$'Pr(>F)'[1]   
  df[2,i] <- t.test(myPCOA$scores[indexWT,i], myPCOA$scores[index10,i])$p.value 
  df[3,i] <- t.test(myPCOA$scores[indexPre,i], myPCOA$scores[indexPost,i])$p.value
}
print(df)

