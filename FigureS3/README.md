---
Title: The role of IAA and characterization of PIN transporters in complex streptophyte alga Chara braunii
Author: "Katarina Kurtović, Stanislav Vosolsobě, Daniel Nedvěd, Karel Müller, Petre Dobrev, Vojtěch Schmidt, Piotr Piszczek, Andre Kuhn, Adrijana Smoljan, Tom J. Fisher, Dolf Weijers, Jiří Friml, John L. Bowman & Jan Petrášek"
output:
  html_document:
  df_print: paged
word_document: default
pdf_document: default
---
  
# AUXIN ACCUMULATION IN *CHARA BRAUNII* UPON EXOGENOUS IAA TREATMENT
  
### R code by Stanislav Vosolsobě

## Loading a libraries

```r
library(emmeans)
```

## Fig. S3 c) Plot of IAA content upon the treatment

```r
metabB <- metab[metab$material=="biomass",][-c(7,35),]
# row 7 removed, to hight
variB <- as.factor(paste(metabB$treatment,metabB$time,sep = " "))
variB <- factor(variB,levels = levels(variB)[c(1,3,2,4,6,5,7,9,8)])
pdf("Metabolism_suppl.pdf",width = 7,height = 7)
par(las=1,mar=c(4,4,2,2))
boxplot(log(metabB$IAA) ~ variB, outline=F,col= c(rep("#efa68eff",3),rep("#41b1a9ff",3),rep("grey",3)),xaxt="n",xlab="Time [h]",ylab="IAA Concentration [pmol/g]")
axis(side=1,1:9,labels = rep(c(1,6,24),3))
legend(x = 1,y = 7,legend = c("10 uM NPA","1 uM IAA","DMSO"),fill = c("#efa68eff","#41b1a9ff","grey"))
points(log(metabB$IAA) ~ jitter(as.numeric(variB),amount = 0.3),pch=16,col="black",cex=1.5)
dev.off()
```

![Boxplot of IAA content](Metabolism_suppl.pdf)

## Statistical analysis of IAA content upon the treatment

```r
m1 <- lm(log(metabB$IAA) ~ variB)
anova(m1)
# p = 2.614e-08 ***

emm <- emmeans(m1, ~ variB)
pairs(emm)

## RESULTS
contrast                 estimate    SE df t.ratio p.value
10uM NPA 1 - 10uM NPA 6   -0.4513 0.581 34  -0.777  0.9968
10uM NPA 1 - 10uM NPA 24   1.9143 0.581 34   3.297  0.0514
#10uM NPA 1 - 1uM IAA 1    -2.9900 0.581 34  -5.150  0.0003
10uM NPA 1 - 1uM IAA 6    -2.4306 0.616 34  -3.947  0.0100
10uM NPA 1 - 1uM IAA 24   -1.1732 0.581 34  -2.021  0.5405
#10uM NPA 1 - DMSO 1        0.3161 0.616 34   0.513  0.9998
10uM NPA 1 - DMSO 6       -0.4329 0.581 34  -0.746  0.9976
10uM NPA 1 - DMSO 24       0.5573 0.581 34   0.960  0.9871
10uM NPA 6 - 10uM NPA 24   2.3656 0.581 34   4.074  0.0071
10uM NPA 6 - 1uM IAA 1    -2.5387 0.581 34  -4.372  0.0031
#10uM NPA 6 - 1uM IAA 6    -1.9792 0.616 34  -3.214  0.0624
10uM NPA 6 - 1uM IAA 24   -0.7219 0.581 34  -1.243  0.9405
10uM NPA 6 - DMSO 1        0.7675 0.616 34   1.246  0.9397
#10uM NPA 6 - DMSO 6        0.0184 0.581 34   0.032  1.0000
10uM NPA 6 - DMSO 24       1.0087 0.581 34   1.737  0.7203
10uM NPA 24 - 1uM IAA 1   -4.9043 0.581 34  -8.447  <.0001
10uM NPA 24 - 1uM IAA 6   -4.3449 0.616 34  -7.055  <.0001
#10uM NPA 24 - 1uM IAA 24  -3.0875 0.581 34  -5.318  0.0002
10uM NPA 24 - DMSO 1      -1.5982 0.616 34  -2.595  0.2259
10uM NPA 24 - DMSO 6      -2.3472 0.581 34  -4.043  0.0077
#10uM NPA 24 - DMSO 24     -1.3570 0.581 34  -2.337  0.3495
1uM IAA 1 - 1uM IAA 6      0.5595 0.616 34   0.908  0.9910
1uM IAA 1 - 1uM IAA 24     1.8168 0.581 34   3.129  0.0757
#1uM IAA 1 - DMSO 1         3.3062 0.616 34   5.368  0.0002
1uM IAA 1 - DMSO 6         2.5571 0.581 34   4.404  0.0029
1uM IAA 1 - DMSO 24        3.5474 0.581 34   6.110  <.0001
1uM IAA 6 - 1uM IAA 24     1.2574 0.616 34   2.042  0.5270
1uM IAA 6 - DMSO 1         2.7467 0.649 34   4.231  0.0046
#1uM IAA 6 - DMSO 6         1.9976 0.616 34   3.244  0.0583
1uM IAA 6 - DMSO 24        2.9879 0.616 34   4.852  0.0008
1uM IAA 24 - DMSO 1        1.4893 0.616 34   2.418  0.3069
1uM IAA 24 - DMSO 6        0.7403 0.581 34   1.275  0.9318
#1uM IAA 24 - DMSO 24       1.7305 0.581 34   2.980  0.1050
DMSO 1 - DMSO 6           -0.7490 0.616 34  -1.216  0.9472
DMSO 1 - DMSO 24           0.2412 0.616 34   0.392  1.0000
DMSO 6 - DMSO 24           0.9902 0.581 34   1.705  0.7391
```


