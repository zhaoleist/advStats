# lab 8
rm(list=ls())

par(mfrow=c(2,2))
# 1
data <- read.delim("clipboard")
myLm <- lm(data[,1] ~ data[,2])
plot(data[,2], data[,1])
abline(myLm$coefficients[1], myLm$coefficients[2], col="red")
anova(myLm)$"Pr(>F)"[1]    # 4.148402e-28, the p-value for the null hypothesis that the two prices are not associated
summary(myLm)$r.squared    # 0.9970831, the r-squared for the association

# 2
path1 <- "http://afodor.github.io/classes/stats2015/caseControlData.txt"
path2 <- "http://afodor.github.io/classes/stats2015/BMI_Data.txt"
dataOTU <- read.table(path1, header=TRUE, sep="\t")
dataBMI <- read.table(path2, header=TRUE, sep="\t", fill=TRUE)    # BMI
# which(count.fields(path2) != 2)
dataBMI <- dataBMI[-(which(is.na(dataBMI[,2]))),]    # remove rows with NA 
  
dataOTU$sample <- gsub("case", "", dataOTU$sample)
dataOTU$sample <- gsub("control", "", dataOTU$sample)
dataOTU$sample <- gsub("_\\d_\\d_\\d*", "", dataOTU$sample)
mergedData <- merge(dataOTU, dataBMI, by.x = "sample", by.y = "studyid")

mergedData_2 <- mergedData[,-1]    # remove first column "sample id" 
pVals <- c()
for (i in 1:(ncol(mergedData_2)-1)){
  myLm <- lm(mergedData_2$bmi ~ mergedData_2[,i])
  pVals[i] <- anova(myLm)$"Pr(>F)"[1]
}
hist(pVals, breaks=20)
# sum(pVals < 0.1)/length(pVals) = 0.1047619
padj <- p.adjust(pVals, method="BH")
if (any(padj < 0.1)==TRUE){
  print (paste("There are", length(which(padj<0.1)), "p-values significant at 10% FDR"))
} else {
  print ("None of these associations significant")
}
# ANSWER:
# The p-values seem to be uniformly distributed.
# The microbial community appears not to be influencing body weight in this cohort.
