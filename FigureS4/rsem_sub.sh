#!/bin/bash
#PBS -N blast_sub
#PBS -l walltime=24:00:00
#PBS -j oe

#######################################################
###### Tento skript stáhne blastuje RNA-Seq db ########
#######################################################

# !!!! Obecně stačí 10 GB RAM a 100 GB HDD
# run: qsub -N mapping -l select=1:ncpus=${thr}:mem=500gb -l walltime=24:00:00 -j oe -v dir=${dir},ncpu=${ncpu} -q ibot star_sub.sh


i=$PBS_ARRAY_INDEX

module add rsem/1.3.3

if [ "${manifest}" = "PE" ]
then
	echo "PE"
	reads=`sed -n "${i}p" ${dir}/manifest${manifest}${suffix}.tsv | cut -f 1,2 --output-delimiter=" "`
    pe="--paired-end"
else
	echo "SI"
	reads=`sed -n "${i}p" ${dir}/manifest${manifest}${suffix}.tsv | cut -f 1`
    pe=""
fi

echo $i
echo $dir
echo $reads
ID=`sed -n "${i}p" ${dir}/manifest${manifest}${suffix}.tsv | cut -f 3`
echo $ID

echo $thr

m1=`grep -v "^N_" ${dir}/star2/${ID}_${manifest}_ReadsPerGene.out.tab | cut -f 3 | paste -s -d "+" | bc`
m2=`grep -v "^N_" ${dir}/star2/${ID}_${manifest}_ReadsPerGene.out.tab | cut -f 4 | paste -s -d "+" | bc`
dd=`echo "${m1}/(${m2}+1)" | bc -l`
if echo "$dd > 10" | bc -l | grep -q 1
then
    strand=1
else
    if echo "$dd < 0.1" | bc -l | grep -q 1
    then
        strand=0
    else
        strand=0.5
    fi
fi
echo $dd $strand
echo $pe

### RSEM quantification on STAR alignment
rsem-calculate-expression --bam --no-bam-output -p ${thr} ${pe} --forward-prob ${strand} ${dir}/star2/${ID}_${manifest}_Aligned.toTranscriptome.out.bam ${dir}/rsemind/rsemind ${dir}/rsem/${ID}_${manifest}


