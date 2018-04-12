# lab 10
rm(list=ls())
library(ggplot2)
# library(gridExtra)
# par(mfrow=c(m,n)) does not work with ggplots. To have multiple graphs on one page or 
# generate multi-page pdf file, package "gridExtra" is needed.

path <- "http://afodor.github.io/classes/stats2015/qPCRWithSampleDays.txt"
myT <- read.table(path, header=TRUE, sep="\t")

# (1)
myT$treatmentStatus <- factor(myT$treatmentStatus)
myT$treatmentStatus <- relevel(myT$treatmentStatus, ref="Treatment")

g <-    # create ggplot object
  ggplot(myT,aes(x=sampleDays, y=Log16S)) + geom_point(aes(col=treatmentStatus), size=2) +
  scale_color_manual(values=c("red","blue","orange", "green"))
g

# (2)

# 8 parameter model
myLm_8p <- lm(myT$Log16S ~ myT$sampleDays * myT$treatmentStatus, x=TRUE, y=TRUE)
myCoefs_8p <- coef(myLm_8p)
g_8p <- g + 
  geom_abline(intercept = myCoefs_8p[1], slope = myCoefs_8p[2], col="red", size=0.6) +
  geom_abline(intercept = myCoefs_8p[1] + myCoefs_8p[3] + myCoefs_8p[6], slope = myCoefs_8p[2] + myCoefs_8p[6], col="blue", size=0.6) +
  geom_abline(intercept = myCoefs_8p[1] + myCoefs_8p[4] + myCoefs_8p[7], slope = myCoefs_8p[2] + myCoefs_8p[7], col="orange", size=0.6) +
  geom_abline(intercept = myCoefs_8p[1] + myCoefs_8p[5] + myCoefs_8p[8], slope = myCoefs_8p[2] + myCoefs_8p[8], col="green", size=0.6) +
  annotate("text", label="8 parameter model", x=750, y=7.5, size=3) 
g_8p

# 5 parameter model:
myLm_5p <- lm(myT$Log16S ~ myT$sampleDays + myT$treatmentStatus, x=TRUE)
myCoefs_5p <- coef(myLm_5p)
g_5p <- g +
  geom_abline(intercept = myCoefs_5p[1], slope = myCoefs_5p[2], col="red", size=0.6) +
  geom_abline(intercept = myCoefs_5p[1] + myCoefs_5p[3], slope = myCoefs_5p[2], col="blue", size=0.6) +
  geom_abline(intercept = myCoefs_5p[1] + myCoefs_5p[4], slope = myCoefs_5p[2], col="orange", size=0.6) +
  geom_abline(intercept = myCoefs_5p[1] + myCoefs_5p[5], slope = myCoefs_5p[2], col="green", size=0.6) +
  annotate("text", label="5 parameter model", x=750, y=7.5, size=3) 
g_5p

# 2 parameter model
myLm_2p <- lm(myT$Log16S ~ myT$sampleDays, x=TRUE)
myCoefs_2p <- coef(myLm_2p)
g_2p <- g +
  geom_abline(intercept = myCoefs_2p[1], slope = myCoefs_2p[2], col="black", size=0.6) +
  annotate("text", label="2 parameter model", x=750, y=7.5, size=3) 
g_2p

# plots <- marrangeGrob(list(g,g_8p,g_5p,g_2p), nrow=1, ncol=1)
# ggsave(plots, filename = "Lab10_Plots_Lei.pdf")

summary(myLm_8p)

# ANSWER:
# The 8 parameter model is the most appropriate. When checking the summary(myLm_8p), it's clear that
# the p-values for three interaction terms are all < 0.05. Thus it's not appropriate to drop these interaction
# terms since they make differences.

# I think the treatment had no statistically significant effect on Log16S. Here we used "Treatment" 
# as background and p-values of "Before Treatment" and "Recovery" are both < 0.05. Before treatment,
# we saw more bugs and when antibiotics were stopped, the bugs went back. Also, from the 8 parameter
# graph, the regression lines of "Before Treatment" and "Recovery" had no clear seperation. From the 5 parameter
# model, the two lines almost completely overlapped. So I would say treatment did not make a difference. 
