# lab 9
rm(list=ls())
myT<-read.table("nc101_scaff_dataCounts.txt",sep="\t",header=TRUE)
myTNorm <- myT
for ( i in 2:ncol(myT))
{
  colSum = sum(myT[,i])
  myTNorm[,i] =myTNorm[,i]/colSum
}

conditions <- factor( c(rep("W02",3),rep("W12",3),rep("W20",5)))

pVals <- vector()
for (i in 1:nrow(myTNorm)){
  myLm <- lm(as.numeric(myTNorm[i,2:12]) ~ conditions, x=TRUE)
  pVals[i] <- anova(myLm)$"Pr(>F)"[1]
}

hist(pVals, breaks=20)    # not uniformaly distributed
print("The p-values are not uniformly distributed.")
pValsAdj <- p.adjust(pVals, method="BH")
sigIndices <- which(pValsAdj < 0.01)    # row numbers of these with significant p-values
print(paste("There are", length(sigIndices), "genes are significant at a BH FDR p-value of 0.01"))

# A new data frame 'myTNorm_sig' only has genes with significant p-values and last column is adjusted p-values
myTNorm_sig <- myTNorm[sigIndices,]
myTNorm_sig$pValues <- pValsAdj[sigIndices]    

# This data frame 'df' has been orderd by significance of p-values
df <- myTNorm_sig[order(myTNorm_sig$pValues),]

############################    ggplot solution:
# library(ggplot2)
# pdf(file="3_Lab9_Plos_Lei.pdf")
# for (i in 1:nrow(df)){
#   values <- as.numeric(df[i,2:12])
#   df_plot <- data.frame(values=values,conditions=conditions)
#   pV <- format(df[i,13], digit=3) 
#   g <- 
#     ggplot(df_plot) + 
#     geom_boxplot(aes(x=conditions, y=values, fill=conditions)) +
#     geom_point(aes(x=conditions, y=values, fill=conditions)) +
#     labs(title=paste(df[i,1], ",", "p = ", pV))
#   print(g)
# }
# dev.off()

############################    Base R solution:
pdf("Lab9_Plots_Lei.pdf")
par(mfrow=c(2,2))
for (i in 1:nrow(df)){
  values <- as.numeric(df[i,2:12])
  pV <- format(df[i,13], digit=3) 
  df_plot <- data.frame(values=values, conditions=conditions)
  plot(df_plot$values ~ df_plot$conditions, col=c("red", "blue", "green"),
       xlab="conditions", ylab=df[i,1], main=paste0(df[i,1], ", p = ", pV))
  stripchart(values ~ conditions, data = df_plot, vertical = TRUE, 
             pch = 21, add=TRUE, col="orange")					
}
dev.off()
