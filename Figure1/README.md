---
Title: The role of IAA and characterization of PIN transporters in complex streptophyte alga Chara braunii
Author: "Katarina Kurtović, Stanislav Vosolsobě, Daniel Nedvěd, Karel Müller, Petre Dobrev, Vojtěch Schmidt, Piotr Piszczek, Andre Kuhn, Adrijana Smoljan, Tom J. Fisher, Dolf Weijers, Jiří Friml, John L. Bowman & Jan Petrášek"
output:
  html_document:
  df_print: paged
word_document: default
pdf_document: default
---
  
# AUXIN EFFECT ON *CHARA BRAUNII* SHOOT REGENERATION
  
### R code by Stanislav Vosolsobě

## Loading a libraries

```r
library(lme4)
library(emmeans)
```

## Fig. 1 a) Analysis of the length of regenerated shoots

```r
trl <- read.table("Treatments_all_length", header = T, sep = "\t")

l1 <- glmer(data = trl, Length ~ Treatment + (1|Replicate), family = "Gamma")
l0 <- glm(data = trl, Length ~ 1 + (1|Replicate), family = "Gamma")

la <- anova(l1,l0, test = "LRT")
# p = 3.753e-12 ***


lm <- emmeans(l1, ~Treatment|Replicate, type = "response")
pairs(lm)
plot(lm)
lci <- confint(lm)
emvar <- lci$Treatment
evars <- unique(trl$Treatment)

pdf("Treatments_all_length.pdf",width=7,height = 7)
par(las=2,mar=c(10,4,4,2))
plot(1,type="n",xlim=c(0,length(emvar)),ylim=c(0,max(trl$Length,na.rm=T)),xaxt="n",xlab="", ylab=expression("Length of new shoots"),main=paste("p =",prettyNum(la$`Pr(>Chisq)`[2],digits=3)))
xp <- (1:length(emvar))-0.5
xo <- match(evars,emvar)

axis(side = 1, at = xp, labels = emvar[xo])
col <- c("gold","skyblue","red","darkgreen")
for(v in 1:length(evars)){
  y=trl$Length[trl$Treatment==evars[v]]
  points(y,x=jitter(rep(xp[v],length(y)),amount = 0.2),pch=16,col="grey")
}
segments(xp,lci$asymp.LCL[xo],xp,lci$asymp.UCL[xo], col=col,lwd=3)
points(lci$response[xo]~xp,pch=23, cex=1.5, col=col, bg= "white",lwd=3)
#segments(xp-0.1,bci$rate[xo],xp+0.1,bci$rate[xo], col=col,lwd=5)
legend("topleft",horiz=F,border = NA, seg.len=4, fill = col, legend = emvar[xo],inset = 0.02,bg = NA,bty = "n")
dev.off()
```

![Boxplot of shoot lengths](Treatments_all_length.pdf)



## Fig. 1 b) Analysis of the number of regenerated shoots

```r
trt <- read.table("Treatments_all_shoots", header = T, sep = "\t")

t1 <- glmer(data = trt, Shoots ~ Treatment + (1|Replicate), family = "poisson")
t0 <- glmer(data = trt, Shoots ~ 1 + (1|Replicate), family = "poisson")

ta <- anova(t1,t0, test = "LRT")
# p = 0.0008021 *** 

tm <- emmeans(t1, ~Treatment, type = "response")
pairs(tm)
plot(tm)
tci <- confint(tm)
emvar <- tci$Treatment
evars <- unique(trt$Treatment)

pdf("Treatments_all_shoots.pdf",width=7,height = 7)
par(las=2,mar=c(10,4,4,2))
plot(1,type="n",xlim=c(0,length(emvar)),ylim=c(0,max(trt$Shoots)),xaxt="n",xlab="", ylab=expression("Number of bottom shoots"),main=paste("p =",prettyNum(ta$`Pr(>Chisq)`[2],digits=3)))
xp <- (1:length(emvar))-0.5
xo <- match(evars,emvar)

axis(side = 1, at = xp, labels = emvar[xo])
col <- c("gold","skyblue","red","darkgreen")
for(v in 1:length(evars)){
  y=trt$Shoots[trt$Treatment==evars[v]]
  points(y,x=jitter(rep(xp[v],length(y)),amount = 0.2),pch=16,col="grey")
}
segments(xp,tci$asymp.LCL[xo],xp,tci$asymp.UCL[xo], col=col,lwd=3)
points(tci$rate[xo]~xp,pch=23, cex=1.5, col=col, bg= "white",lwd=3)
#segments(xp-0.1,bci$rate[xo],xp+0.1,bci$rate[xo], col=col,lwd=5)
legend("topleft",horiz=F,border = NA, seg.len=4, fill = col, legend = emvar[xo],inset = 0.02,bg = NA,bty = "n")
dev.off()
```

![Boxplot of shoot numbers](Treatments_all_shoots.pdf)


