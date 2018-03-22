# lab 7
rm(list=ls())
library(DESeq)

# 2
numRows = 3000
numCols = 10
for( i in 1:numCols)
  
  myFrame <- data.frame(1:numRows)

#initiate the data.frame with the correct # of rows to suppress error messages.
#likely, there are much better ways to do this!
names(myFrame)[1] <- "tempColumn"

for( i in 1: numCols)
{
  vals <- vector(length=numRows)
  
  for( j in 1:numRows)
  {
    aMean = j /10
    aMean = max( aMean,5)
    aVar = aMean+ 5* aMean 
    aVal = round( max( rnorm(1,mean=aMean,sd=sqrt(aVar)), 1))
    vals[j] = aVal
  }
  
  colName <- paste( "sample" , i ,sep="")
  
  myFrame[[colName]] = vals
}

myFrame["tempColumn"] <- NULL
row.names(myFrame) <- paste("Gene_",1:numRows,sep="")

pdf("Lab7_Plots_Lei.pdf")
par(mfrow=c(2,2))
# 3
cds <- newCountDataSet(myFrame, c(rep("A", 5), rep("B", 5)))
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds_2 <- nbinomTest(cds,"A","B") 

means <- apply(counts(cds,normalized=TRUE), 1,mean)
myInfo <- fitInfo(cds)

# DeSeq's estimate of dispersion calculated for each gene
plot(means, means * means* myInfo$perGeneDispEsts)

# the fit that DeSeq uses for the pooled "dispersion"
lines(means, means*means* myInfo$dispFunc(means),col="RED")

# the dispersion for each gene that DeSeq actually uses for inference
points(means, means * means* fData(cds)[,1], col="YELLOW")
hist(cds_2$pval, breaks = 20, main="Histogram of Question (3)")
sum(cds_2$padj < 0.1)

# ANSWER (3):
# The variance DeSeq uses for inference (yellow points) are the maximun values of empirical
# values and fitted values. That's why they are "above" the fitted line. By this way, it avoids
# some situations that the dispersion of a gene is underestimated. Thus, it's conservative.

# Basically it's not that close to uniform with fewer p-values less than 0.1 and more p-values close to 1.
# Yes, this data analysis path through DeSeq is slightly conservative and I found zero gene 
# significantly differernt at 10% FDR.

# 4
cds <- newCountDataSet(myFrame, c(rep("A", 5), rep("B", 5)))
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds, sharingMode="gene-est-only")
cds_3 <- nbinomTest(cds,"A","B")

means <- apply(counts(cds,normalized=TRUE), 1,mean)
myInfo <- fitInfo(cds)

# DeSeq's estimate of dispersion calculated for each gene
plot(means, means * means* myInfo$perGeneDispEsts)

# the fit that DeSeq uses for the pooled "dispersion"
lines(means, means*means* myInfo$dispFunc(means),col="RED")

# the dispersion for each gene that DeSeq actually uses for inference
points(means, means * means* fData(cds)[,1], col="YELLOW")
hist(cds_3$pval, breaks=20, main="Histogram of Question (4)")
sum(cds_3$padj < 0.1)

# ANSWER (4):
# Now the two (yellow and black points) are highly similar. The p-values are still not 
# close to uniform with more p-values less than 0.05. 
# I see 28 genes are significantly different at 10% FDR . 
# This is less conservative than the default option.  

# 5
pVals <- apply(myFrame, 1, function(x) {t.test(x[1:5], x[6:10], conf.level = 0.9)$p.value})
# pVals <-vector()
# for (i in 1:nrow(myFrame)){
#   pVals[i]  <- t.test(myFrame[i,1:5], myFrame[i, 6:10])$p.value
# }
hist(pVals, breaks=20, main="Histogram of p-values from t-test")
plot(pVals, cds_2$pval, main="p-values: t-test vs. DeSeq(default)")
plot(pVals, cds_3$pval, main="p-values: t-test vs. DeSeq(gene-est-only)")

# ANSWER (5)
# The most conservative test is the DeSeq with the default choice.
# T-test comes cloest to a uniform distribution of p-values.
dev.off()

