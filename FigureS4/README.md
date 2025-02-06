---
Title: The role of IAA and characterization of PIN transporters in complex streptophyte alga Chara braunii
Author: "Katarina Kurtović, Stanislav Vosolsobě, Daniel Nedvěd, Karel Müller, Petre Dobrev, Vojtěch Schmidt, Piotr Piszczek, Andre Kuhn, Adrijana Smoljan, Tom J. Fisher, Dolf Weijers, Jiří Friml, John L. Bowman & Jan Petrášek"
output:
  html_document:
  df_print: paged
word_document: default
pdf_document: default
---
  
# ANALYSIS OF RNA-seq LIBRARIES



## Import of dataset and its adaptation

```r
profile <- read.table("profile",header=T)
db.list <- read.table("RNAseq.db.list",header=T)

# Gene definition
PINs <- c("PINa","PINb","PINc","PINd","PINe","PINf")
PINacc <- c("GBG79698","GBG79697","GBG64280","GBG90244","GBG90247","GBG42734")
PINgid <- c("CBR_g29962","CBR_g29961","CBR_g41200","CBR_g50423","CBR_g50425","CBR_g84230")
PINacc.1 <- paste(PINacc,".1", sep="")
PINcols <- c("red","orange","green3","blue3","purple","brown")

# New variables
acc <- profile$gene_id
dbs <- unlist(lapply(strsplit(colnames(profile)[-1],split="\\."),"[",1))
dbs_id <- match(dbs,db.list$Accession)
vars <- db.list$Sample[dbs_id]
```

## Barplot of expression, filtered TPM > 0.5

```r
pdf("profiles_veg_ncbi_0.5.pdf",width=12,height=6)
par(las=2,mar=c(6,4,3,1))
dbs_id <- match(dbs,db.list$Accession)
vars <- db.list$Sample[dbs_id]
sel <- (db.list$Strain[dbs_id]=="S276")&(db.list$Treatment[dbs_id]=="WT")&(db.list$Tissue[dbs_id]!="nodal")&(db.list$Source[dbs_id]=="NCBI")
vars2 <- db.list$Tissue[dbs_id]
setdiff(colnames(profile[,-1][,sel]),db.list$Accession[dbs_id][sel]) # are identical
plot(1,type="n",xlim=c(1-0.3,length(unique(vars2[sel]))+0.3), ylim=c(0,sqrt(20)), xaxt="n",xlab="",ylab="TPM",yaxt="n")
axis(side=1,at=unique(as.numeric(as.factor(vars2[sel]))),labels=unique(vars2[sel]))
axis(side=2,at=(0:5),labels=(0:5)^2)
abline(v=0:length(unique(vars2))+0.5,col="grey")
for(i in 1:length(PINgid)){
  k <- which(acc==PINgid[i])
  print(k)
  y <- profile[k,-1][sel]
  y[y<0.5] <- 0
  rect(ytop=sqrt(tapply(unlist(y),vars2[sel],median)), ybottom=0, xleft=(1:length(unique(vars2[sel])))-0.35+0.1*i-0.04,  xright=(1:length(unique(vars2[sel])))-0.35+0.1*i+0.04, border=NA, col=alpha(PINcols[i],alpha=0.5))
  points(sqrt(unlist(y))~I(jitter(as.numeric(as.factor(vars2[sel])),factor=0.1)-0.35+0.1*i),pch=16,col=PINcols[i])
  #segments(y0=sqrt(tapply(unlist(profile[k,-1][sel]),varss[sel],median)),  y1=sqrt(tapply(unlist(profile[k,-1][sel]),varss[sel],median)),  x0=(1:length(varss))-0.35+0.1*i-0.06,  x1=(1:length(varss))-0.35+0.1*i+0.06, lwd=3, col=PINcols[i])
}
#legend(x=0.5,y=sqrt(9),fill=PINcols,legend=paste(PINs,PINacc),horiz=T,border=F,bg="white")
par(xpd=T)
legend("topleft",inset = c(0,-0.12),fill=PINcols,legend=paste(PINs),horiz=T,border=F,bg="white")
par(xpd=F)
dev.off()
```





