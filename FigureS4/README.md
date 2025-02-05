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

Quantification of RNA-seq libraries were computed on METACENRUM computer cluster. BASH scripts are presented.

### Path specification

```sh
dir="/storage/brno3-cerit/home/vosolsob/Exprese"
cd $dir
```

## Mapping of RNA-seq reads by STAR to genome of *Chara braunii*

### genome index for STAR

```sh
qsub -N gind -l select=1:ncpus=40:mem=50gb -l walltime=24:00:00 -j oe -v dir=${dir} -q ibot
echo $dir
module add star-2.7.7a
STAR --runThreadN 40 --runMode genomeGenerate --genomeDir ${dir}/genind --genomeFastaFiles ${dir}/genome/GCA_003427395.1_Cbr_1.0_genomic.fna --sjdbGTFfile ${dir}/genome/GCA_003427395.1_Cbr_1.0_genomic.gtf
```

### Alignment of reads by STAR

```sh
thr=10

# first dataset
gz="no"
suffix=""

manifest="PE"
ndb=`grep -c . manifest${manifest}.tsv`
qsub -J 1-${ndb} -N mapping -l select=1:ncpus=${thr}:mem=50gb -l walltime=24:00:00 -j oe -v dir=${dir},thr=${thr},manifest=${manifest},gz=${gz},suffix=${suffix} star_sub-2.0.sh
manifest="SI"
ndb=`grep -c . manifest${manifest}.tsv`
qsub -J 1-${ndb} -N mapping -l select=1:ncpus=${thr}:mem=50gb -l walltime=24:00:00 -j oe -v dir=${dir},thr=${thr},manifest=${manifest},gz=${gz},suffix=${suffix} star_sub-2.0.sh
# second dataset
gz="yes"
suffix="_orig"

manifest="PE"
ndb=`grep -c . manifest${manifest}${suffix}.tsv`
qsub -J 1-${ndb} -N mapping -l select=1:ncpus=${thr}:mem=50gb -l walltime=24:00:00 -j oe -v dir=${dir},thr=${thr},manifest=${manifest},gz=${gz},suffix=${suffix} star_sub-2.0.sh
manifest="SI"
ndb=`grep -c . manifest${manifest}${suffix}.tsv`
qsub -J 1-${ndb} -N mapping -l select=1:ncpus=${thr}:mem=50gb -l walltime=24:00:00 -j oe -v dir=${dir},thr=${thr},manifest=${manifest},gz=${gz},suffix=${suffix} star_sub-2.0.sh```
```

## Quantification of STAR alignments by RSEM

### genome index for RSEM

```sh
qsub -N trind
module add rsem/1.3.3
rsem-prepare-reference --gtf ${dir}/genome/GCA_003427395.1_Cbr_1.0_genomic.gtf ${dir}/genome/GCA_003427395.1_Cbr_1.0_genomic.fna ${dir}/rsemind/rsemind
```

### Quantification by RSEM

```sh
thr=4

# first dataset
gz="no"
suffix=""

manifest="PE"
ndb=`grep -c . manifest${manifest}.tsv`
qsub -J 1-${ndb} -N quant -l select=1:ncpus=${thr}:mem=50gb -l walltime=24:00:00 -j oe -v dir=${dir},thr=${thr},manifest=${manifest},gz=${gz},suffix=${suffix} rsem_sub.sh

manifest="SI"
ndb=`grep -c . manifest${manifest}.tsv`
qsub -J 1-${ndb} -N quant -l select=1:ncpus=${thr}:mem=50gb -l walltime=24:00:00 -j oe -v dir=${dir},thr=${thr},manifest=${manifest},gz=${gz},suffix=${suffix} rsem_sub.sh

# second dataset
gz="yes"
suffix="_orig"

manifest="PE"
ndb=`grep -c . manifest${manifest}${suffix}.tsv`
qsub -J 1-${ndb} -N quant -l select=1:ncpus=${thr}:mem=50gb -l walltime=24:00:00 -j oe -v dir=${dir},thr=${thr},manifest=${manifest},gz=${gz},suffix=${suffix} rsem_sub.sh

manifest="SI"
ndb=`grep -c . manifest${manifest}${suffix}.tsv`
qsub -J 1-${ndb} -N quant -l select=1:ncpus=${thr}:mem=50gb -l walltime=24:00:00 -j oe -v dir=${dir},thr=${thr},manifest=${manifest},gz=${gz},suffix=${suffix} rsem_sub.sh
```

### Merging of individual abundancies quantified by RSEM

```sh
dir="/storage/brno12-cerit/home/vosolsob/Exprese"
cd $dir/rsem

lf=(`ls *genes.results`)
cut -f 1 $lf > profile
for f in ${lf[*]}
do
echo "Processing $f" # always double quote "$f" filename
#wc -l $f
name=`echo $f | cut -d "_" -f 1`
cut -f 6 $f | sed -e "s/TPM/${name}/" | paste profile - > profile2
rm profile
mv profile2 profile
done
```

## Visualisation of transcriptomic profiles in R

### Import of dataset and its adaptation

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

### Barplot of expression, filtered TPM > 0.5

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

![Barplot of PIN gene expression](profiles_veg_ncbi_0.5.pdf)





