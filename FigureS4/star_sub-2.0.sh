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
 
#module add star-2.7.7a
module add star/2.7.10b-gcc-10.2.1-6cnaeq7
#module add samtools/samtools-1.11-intel-19.0.4-wzth4e4
#module add samtools/1.14-intel-19.0.4-nlgg3wv
module add salmon/1.10.1-gcc-10.2.1-5egjrbv

if [ "${manifest}" = "PE" ]
then
	echo "PE"
	reads=`sed -n "${i}p" ${dir}/manifest${manifest}${suffix}.tsv | cut -f 1,2 --output-delimiter=" "`
else
	echo "SI"
	reads=`sed -n "${i}p" ${dir}/manifest${manifest}${suffix}.tsv | cut -f 1`
fi

echo $i
echo $dir
echo $reads
ID=`sed -n "${i}p" ${dir}/manifest${manifest}${suffix}.tsv | cut -f 3`
echo $ID

echo $thr


if [ "${gz}" = "yes" ]
then
	echo "zipped"
    read="--readFilesCommand gunzip -c"
else
	echo "unzipped"
	read=""
fi
echo $read

#--sjdbOverhang 149
STAR --runThreadN ${thr} --twopassMode Basic --sjdbScore 2 --limitSjdbInsertNsj 1000000 --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 10000000 --alignMatesGapMax 1000000 --outSAMunmapped Within --outFilterType BySJout --outSAMattributes NH HI AS NM MD --quantMode TranscriptomeSAM GeneCounts --sjdbGTFfile ${dir}/genome/GCA_003427395.1_Cbr_1.0_genomic.gtf --genomeDir ${dir}/genind --readFilesIn ${reads} --outFileNamePrefix ${dir}/star2/${ID}_${manifest}_ -outSAMtype BAM SortedByCoordinate ${read}


### SALMON quantification on STAR alignment

salmon quant --targets ${dir}/genome/GCA_003427395.1_Cbr_1.0_cds_from_genomic_short_names.fna --libType ISR --gencode --gcBias --seqBias --posBias -p ${thr} --numBootstraps 30 -o ${dir}/salmon/${ID}_${manifest}.tpm --alignments ${dir}/star2/${ID}_${manifest}_Aligned.toTranscriptome.out.bam -g ${dir}/genome/GCA_003427395.1_Cbr_1.0_genomic.gtf



