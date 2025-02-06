---
Title: The role of IAA and characterization of PIN transporters in complex streptophyte alga Chara braunii
Author: "Katarina Kurtović, Stanislav Vosolsobě, Daniel Nedvěd, Karel Müller, Petre Dobrev, Vojtěch Schmidt, Piotr Piszczek, Andre Kuhn, Adrijana Smoljan, Tom J. Fisher, Dolf Weijers, Jiří Friml, John L. Bowman & Jan Petrášek"
output:
  html_document:
  df_print: paged
word_document: default
pdf_document: default
---
  
# COMPLEMENTATION OF *pin2* MUTATION BY CbPINa AND CbPINc IN *ARABIDOPSIS THALIANA*

## Loading a libraries

```r
require(CircStats)
require(circular)
require(bpnreg)
require(plotrix)
require(boot)
```

## Loading a data

```r
comp <- read.table("complementation", header=T)
```

## Function, which compute circular stats within the bootstrapping procedure

```r
circ_boot <- function(dta,i){
  m <- tapply(dta$angle[i],dta$genotype[i],mean.circular)
  r <- tapply(dta$angle[i],dta$genotype[i],rho.circular)
  return(c(m,r))
}
```

## Boostrapping and estimation of CIs

```r
bb <- boot(dta,circ_boot,R=1000,strata = dta$genotype,stype = "i")

means.ci <- matrix(NA,nrow = nlevels(dta$genotype)+1,ncol=3)
rownames(means.ci) <- c(names(bb$t0)[1:nlevels(dta$genotype)],"random")
colnames(means.ci) <- c("LL","AVG","UL")
for(i in 1:nlevels(dta$genotype)){
  means.ci[i,c(1,3)] <- quantile.circular(circular(bb$t[,i],type = "angles",units = "degrees",rotation = "clock",zero = -pi/2),probs = c(0.025,0.975))
  means.ci[i,2] <- mean.circular(circular(bb$t[,i],type = "angles",units = "degrees",rotation = "clock",zero = -pi/2))
}

rho.ci <- matrix(NA,nrow = nlevels(dta$genotype)+1,ncol=3)
rownames(rho.ci) <- c(names(bb$t0)[nlevels(dta$genotype)+1:nlevels(dta$genotype)],"random")
colnames(rho.ci) <- c("LL","AVG","UL")
for(i in nlevels(dta$genotype)+1:nlevels(dta$genotype)){
  rho.ci[i-nlevels(dta$genotype),c(1,3)] <- quantile(bb$t[,i],probs = c(0.025,0.975))
  rho.ci[i-nlevels(dta$genotype),2] <- mean(bb$t[,i])
}
```

## Generate set of random angles and CIs for comparison

```r
randstat <- matrix(NA,ncol=2,nrow=1000)
colnames(randstat) <- c("mean","rho")
for(i in 1:1000){
  samp <- runif(30,min = 0,max = 360)
  randstat[i,1] <- mean.circular(circular(samp,type = "angles",units = "degrees",rotation = "clock",zero = -pi/2))
  randstat[i,2] <- rho.circular(circular(samp,type = "angles",units = "degrees",rotation = "clock",zero = -pi/2))
}

means.ci[nrow(means.ci),c(1,3)] <- quantile.circular(circular(randstat[,1],type = "angles",units = "degrees",rotation = "clock",zero = -pi/2),probs = c(0.025,0.975))
means.ci[nrow(means.ci),2] <- mean.circular(circular(randstat[,1],type = "angles",units = "degrees",rotation = "clock",zero = -pi/2))

rho.ci[nrow(rho.ci),c(1,3)] <- quantile(randstat[,2],probs = c(0.025,0.975))
rho.ci[nrow(rho.ci),2] <- mean(randstat[,2])
```

## Final figure

```r
pdf("complementation_stat.pdf",width = 8, height = 8)
cols <- c("#efa68eff","#41b1a9ff","darkseagreen3","gray25","grey")
par(mar=c(1,1,1,1))
plot.circular(dta$angle,stack = F,cex = 0)
segments(-1,0,1,0)
segments(0,-1,0,1)
for(t in 1:(nrow(means.ci)-1)){
  points.circular(dta$angle[dta$genotype==rownames(means.ci)[t]],stack = F,col=cols[t],start.sep = t/50)
}
for(r in 1:nrow(means.ci)){
  points.circular(circular(means.ci[r,2],type = "angles",units = "degrees",rotation = "clock",zero = -pi/2),col=cols[r],start.sep = rho.ci[r,2]-1,cex=2)
  lines.circular(x=circular(means.ci[r,c(2,2)],type = "angles",units = "degrees",rotation = "clock",zero = -pi/2),col=cols[r],y=rho.ci[r,c(1,3)]-1,lwd=3)
}
draw.arc(0,0,rho.ci[-nrow(means.ci),2],deg1 = -means.ci[-nrow(means.ci),1]-90, deg2 = -means.ci[-nrow(means.ci),3]-90,col=cols[-nrow(means.ci)],lwd=3)
draw.arc(0,0,rho.ci[nrow(means.ci),2],deg1 = 0, deg2 = 360,col=cols[nrow(means.ci)],lwd=3)
legend(-0.3,0.8,legend = rownames(means.ci),fill = cols)
dev.off()
```

![Circular diagram](complementation_stat.pdf)








