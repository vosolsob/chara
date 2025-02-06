---
Title: The role of IAA and characterization of PIN transporters in complex streptophyte alga Chara braunii
Author: "Katarina Kurtović, Stanislav Vosolsobě, Daniel Nedvěd, Karel Müller, Petre Dobrev, Vojtěch Schmidt, Piotr Piszczek, Andre Kuhn, Adrijana Smoljan, Tom J. Fisher, Dolf Weijers, Jiří Friml, John L. Bowman & Jan Petrášek"
output:
  html_document:
  df_print: paged
word_document: default
pdf_document: default
---
  
# ANALYSIS OF PIN FROM *CHARA BRAUNII* IN TOBACCO BY-2 CELLS
  
### R code by Stanislav Vosolsobě

## Loading a libraries

```r
library(peripheral)   # installed by http://kfrserver.natur.cuni.cz/lide/vosolsob/Peripheral/

library(nlme)
library(lme4)
library(emmeans)
library(scales)

```



## Fig. 5 b) Analysis of PM localisation of CbPINa and CbPINc

Detailed protocol was published in Vosolsobě, S., Schwarzerová, K., Petrášek, J. (2018): Determination of Plasma Membrane Partitioning for Peripherally-associated Proteins. JoVE (Journal of Visualized Experiments) (136), e57837

```r
profile <- read.profile(file="Results.csv")

# Now, the model of PM and cytoplasmic distribution based on previuos analysis is loaded
# For details, follow the protocol in http://kfrserver.natur.cuni.cz/lide/vosolsob/Peripheral/
load("model.RData")

filtered <- profile
md <- distrib(filtered, mp)
mdf <- filter(md, cut.off = 0.2)

df <- data.frame(rownames(mdf$affinity), mdf$affinity)
colnames(df) <- c(" ", colnames(mdf$affinity))
write.csv(df, file = "distrib.csv", row.names = FALSE, quote = FALSE)

pdf("distrib.pdf", width=10, height = 8)
boxplot(mdf, variants=c("var2"))
dev.off()

anova(lm(log(mdf$affinity$PM.aff)~mdf$affinity$var2))
# p = 3.898e-13 ***
```

![Boxplot of PM affinity](distrib.pdf)


## Fig. 5 c) Accumulation of radiolabelled IAA in BY-2 cells expressing CbPINa and CbPINc

### Data preparation

```r
assay0 <- read.table("assay", header = T, stringsAsFactors = T)

# Creation of factor variables
levels(assay0$EST) <- c("Ctr","Ind")
levels(assay0$NPA) <- c("DMSO","NPA")
Repl <- as.factor(paste(assay0$Date,assay0$EST,assay0$NPA,assay0$Replication,sep="-"))
assay <- cbind(assay0,Repl)
assay$pmolPM[160] <- NA
```

# Fit of asymptotic growth curve to the data 

```r
outable <- data.frame(Date=NULL,Day=NULL,Gene=NULL,EST=NULL,NPA=NULL,Asym=NULL,R0=NULL,lrc=NULL)
i <- 1
for(dt in levels(assay$Date)){
  print(dt)
  for(e in levels(assay$EST)){
    for(n in levels(assay$NPA)){
      for(r in unique(assay$Replication)){
      #gassay <- groupedData(pmolPM~Time|Replication,assay[(assay$EST==e)&(assay$NPA==n)&(assay$Date==dt),])
      gassay <- assay[(assay$EST==e)&(assay$NPA==n)&(assay$Date==dt)&(assay$Replication==r),]
      fm0 <- tryCatch(nls(data=gassay,pmolPM ~ SSasymp(Time, Asym, R0, lrc)),error=function(e) NA)
      if(is.na(fm0)){
        print(paste("ERROR in NLS",dt,e,n,r))
        lm0 <- lm(data=gassay,pmolPM ~ Time)
        outable[i,c("Asym","R0","lrc" )] <- c(NA,coef(lm0)[1],log(coef(lm0)[1]))
      }else{
        outable[i,names(coef(fm0))] <- coef(fm0)
      }
      outable[i,"Date"] <- dt
      outable[i,"Day"] <- gassay$Day[1]
      outable[i,"Gene"] <- gassay$Gene[1]
      outable[i,"EST"] <- e
      outable[i,"NPA"] <- n
      i <- i + 1
      }
    }
  }
}
```

### Plot of fitted asymptotic parameter with statistics

```r
pdf("Accumulation.pdf",width=8,height = 7)
  par(mfrow=c(1,1), las=2,mar=c(6,4,4,1))
  d <- 5
  plot(1, type="n",xlim=c(0.5,8.5),ylim=c(0.5,1.7),xaxt="n",xlab="",ylab="")
  k <- 0
  for(g in unique(outable$Gene)){
    subtable <- outable[(outable$Day==d)&(outable$Gene==g),]
    if(length(unique(subtable$Date))>1){
      sm2 <- lmer(data=subtable, Asym~EST*NPA + (1|Date))
      sm1 <- lmer(data=subtable, Asym~EST+NPA + (1|Date))
      sm0 <- lmer(data=subtable, Asym~NPA + (1|Date))
      sta <- paste("EST:NPA = ",prettyNum(anova(sm1,sm0,test="Chisq")$`Pr(>Chisq)`[2],digits=3),
                   ", EST = ",prettyNum(anova(sm2,sm1,test="Chisq")$`Pr(>Chisq)`[2],digits=3), sep= "")
    }else{
      sm2 <- lm(data=subtable, Asym~EST*NPA)
      sm1 <- lm(data=subtable, Asym~EST+NPA)
      sm0 <- lm(data=subtable, Asym~NPA)
      sta <- paste("EST:NPA = ",prettyNum(anova(sm1,sm0,test="Chisq")$`Pr(>Chi)`[2],digits=3),
                   ", EST = ",prettyNum(anova(sm2,sm1,test="Chisq")$`Pr(>Chi)`[2],digits=3), sep= "")
    }
    tm <- emmeans(sm2, ~EST|NPA, type = "response")
    pp <- summary(pairs(tm))$p.value
    pa <- symnum(pp,cutpoints = c(0,  .001,.01,.05, .1, 1),
                 symbols = c(" ***"," **"," *"," .",""))
    sig <- paste(prettyNum(pp,digits=3),pa,sep="")
    #plot(tm,comparisons = T)
    #tci <- confint(tm)
    #emvar <- tci$Treatment
    #evars <- unique(trt$Treatment)
    
    boxplot(data=subtable, Asym~EST*NPA, main=paste(g,"day",d,sta),col=c("#efa68eff","#41b1a9ff"),add=T,at=(1:4)+k)
    # vertical=T,method="jitter",pch=19,col=alpha("black",alpha=0.2))
    dmso.max <- max(subtable[subtable$NPA=="DMSO",]$Asym,na.rm=T)
    npa.max <- max(subtable[subtable$NPA=="NPA",]$Asym,na.rm=T)
    text(x=1.5+k,y=0.15+dmso.max,sig[1])
    segments(x0 = 1+k,y0 = 0.1+dmso.max,x1 = 2+k,y1 = 0.1+dmso.max)
    segments(x0 = 1+k,y0 = 0.1+dmso.max,x1 = 1+k,y1 = 0.05+dmso.max)
    segments(x0 = 2+k,y0 = 0.1+dmso.max,x1 = 2+k,y1 = 0.05+dmso.max)
    
    text(x=3.5+k,y=0.15+npa.max,sig[2])
    segments(x0 = 3+k,y0 = 0.1+npa.max,x1 = 4+k,y1 = 0.1+npa.max)
    segments(x0 = 3+k,y0 = 0.1+npa.max,x1 = 3+k,y1 = 0.05+npa.max)
    segments(x0 = 4+k,y0 = 0.1+npa.max,x1 = 4+k,y1 = 0.05+npa.max)
    k <- k + 4
  }
dev.off()
```

![Boxplot of IAA accumulation](Accumulation.pdf)

