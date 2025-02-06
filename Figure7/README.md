---
Title: The role of IAA and characterization of PIN transporters in complex streptophyte alga Chara braunii
Author: "Katarina Kurtović, Stanislav Vosolsobě, Daniel Nedvěd, Karel Müller, Petre Dobrev, Vojtěch Schmidt, Piotr Piszczek, Andre Kuhn, Adrijana Smoljan, Tom J. Fisher, Dolf Weijers, Jiří Friml, John L. Bowman & Jan Petrášek"
output:
  html_document:
  df_print: paged
word_document: default
pdf_document: default
---
  
# EFFECT OF IAA AND BA TREATMENT OF CYTOPLASMIC STREAMING IN *CHARA BRAUNII*

## Loading a libraries

```r
require(CircStats)
require(circular)
require(bpnreg)
require(plotrix)
require(boot)
```


## Import of dataset and its adaptation

```r
# set calibration (pixels/um)
pxpum <- 1.476 

track <- read.table("BA4",header=T)

# New variables from description of variants
conc <- unlist(lapply(strsplit(track$video_id,split = "_"),`[`,1))
movi <- as.numeric(unlist(lapply(strsplit(track$video_id,split = "_"),`[`,2)))

# Merging of ALL DMSO variant together
conc2 <- conc
conc2[track$treatment=="dmso"] <- "0"
treatment2 <- paste(track$treatment,conc2)
```

## Computing of velocities from original particle coordinates

```r
ps <- track$slice[1]
px <- track$x[1]
py <- track$y[1]
id <- 0
trt <- NULL #treatment
tr2 <- NULL #treatment with merged DMSO
con <- NULL #concentration
mid <- NULL #movie ID
vel <- NULL #velocity
slc <- NULL #slice
pid <- NULL #particle
cel <- NULL #cell ID
j <- 1
for(i in 2:nrow(track)){
  if( (track$slice[i]<ps) | (track$slice[i] - ps > 3)){
    id <- id+1
    px <- track$x[i]
    py <- track$y[i]
    ps <- track$slice[i]
    print("****")
    next
  }
  trt[j] <- track$treatment[i]
  tr2[j] <- treatment2[i]
  con[j] <- conc[i]
  mid[j] <- movi[i]
  vel[j] <- sqrt( (px - track$x[i])^2 + (py - track$y[i])^2 )
  pid[j] <- id
  slc[j] <- track$slice[i]
  cel[j] <- track$cell[i]
  px <- track$x[i]
  py <- track$y[i]
  ps <- track$slice[i]
  if(vel[j]>200){
    print(track[i-1,])
    print(track[i,])
    print(track[i+1,])
    next   #skip this measurement, it will be raplaced by the next one
  }
  j <- j+1
}
trt <- as.factor(trt)
tr2 <- as.factor(tr2)
con <- as.factor(con)
pid <- as.factor(pid)
```

# Computing of mean velocity for each particle

```r
trt3 <- (NULL) #treatment
tr23 <- (NULL) #treatment
con3 <- (NULL) #concentration
mid3 <- (NULL) #movie ID
vel3 <- (NULL) #velocity
slc3 <- (NULL) #slice
pid3 <- (NULL) #particle
cel3 <- (NULL) #cell
k<-1
for(t in 1:nlevels(trt)){
  ti <- which(trt==levels(trt)[t])
  for(c in 1:nlevels(con)){
    ci <-which(con[ti]==levels(con)[c])
    for(p in 1:nlevels(droplevels(pid[ti][ci]))){
      pi <- which(pid[ti][ci]==levels(droplevels(pid[ti][ci]))[p])
      vm <- median(vel[ti][ci][pi])
      trt3[k] <- as.character(trt[ti][ci][1])
      tr23[k] <- as.character(tr2[ti][ci][1])
      con3[k] <- as.character(con[ti][ci][1])
      mid3[k] <- mid[ti][ci][pi][1]
      vel3[k] <- vm
      slc3[k] <- mean(slc[ti][ci][pi])
      cel3[k] <- cel[ti][ci][pi][1]
      pid3[k] <- p
      k <- k+1
    }
  }
}
```

## Final plotting

```r
pdf("End_tracking.pdf",width=6,height=6)
par(mar=c(5,4,1,1),las=2)
boxplot((vel3/pxpum*2)~tr23,col=c("yellowgreen","seagreen","grey","gold","orange"),xlab = "")
dev.off()
```


![Boxplot of streaming](End_tracking.pdf)
