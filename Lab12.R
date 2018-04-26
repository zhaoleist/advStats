# lab 12
rm(list=ls())
library(nlme)

path <- "http://afodor.github.io/classes/stats2015/prePostPhylum.txt"
myT <- read.table(path, header=TRUE, sep="\t")
myT <- myT[myT$time == "POST",]

myColors <- c("darkred", "red1","darkblue", "dodgerblue2", "darkgreen", "green",
              "orange1", "yellow", "deeppink", "pink1", "ivory")
palette(myColors)

coefGLS <- vector()
pValuesMixed <- vector()
index <- 1

pdf("Lab12_Plots_Lei.pdf")
par(mfrow=c(3,2))
for (i in 5:ncol(myT)){
  bug <- myT[,i]
  cage <- myT$cage
  genotype <- myT$genotype
  df <- data.frame(bug, cage, genotype)
    
  M.gls <- gls( bug ~ genotype , method = "REML", correlation = corCompSymm( form = ~ 1 | cage),data=df)
  coefGLS[index] <- coef(M.gls$modelStruct[1]$corStruct,unconstrained=FALSE)[[1]]
    
  M.mixed <- lme( bug ~ genotype, method= "REML", random = ~1 | cage, data = df)
  pValuesMixed[index] <- unclass(summary(M.mixed))$tTable[2,5]
  
  boxplot(myT[,i] ~ myT$cage, boxfill=myColors,
          main=paste0(names(myT)[i], ", rho = ", format(coefGLS[index], digits=4)))
  
  # text(df$cage, par("usr")[3], labels = levels(unique(df$cage)), srt = 45,  adj = c(1.1,1.1), xpd = TRUE, cex=1.2)
  stripchart(bug ~ cage, data = df,vertical = TRUE, pch = 21, add=TRUE)
    
  index <- index + 1
}
dev.off()

#ANWSERS:
# (1) 
# Visually, there seems to be cage effect.
# (2)
if (any(p.adjust(pValuesMixed, method = "BH") < 0.1)){
  indexSig <- which(p.adjust(pValuesMixed, method = "BH") < 0.1)
  print(paste0(names(myT)[4 + indexSig], " is significantly different for genotype in the mixed model at a 10% FDR"))
}