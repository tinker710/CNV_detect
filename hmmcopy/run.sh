#!/bin/bash

#slurm options
#SBATCH -p intel-e5,amd-ep2
#SBATCH -q normal
#SBATCH -J bowtie2_index

#SBATCH -c 6
#SBATCH -o %j.log
#SBATCH --mem 20g

module load bowtie/2.4.2
module load samtools/1.13
module load bwa

sample_name="A81L_FDSW202625676-1r_BHKHF5DSXY-new_L1"
ref_genome=./genome/genome.fa
## 1. build bwa and bowtie2 index if the index files don't exist
#bwa index ./genome/genome.fa
#bowtie2-build ./genome/genome.fa ./genome/genome.fa

# 2. alignment
bwa mem -t 5 /storage/publicdata/ref/bwa/hg38_Ensembl/genome.fa raw_data/${sample_name}_1.fq.gz raw_data/${sample_name}_2.fq.gz | samtools view -Sb > ${sample_name}.bam
samtools sort -@ 4 -m 2G -O bam -o ${sample_name}_sorted.bam A81L_FDSW202625676-1r_BHKHF5DSXY-new_L1.bam

samtools index ${sample_name}_sorted.bam ${sample_name}_sorted.bam.bai

# 3 hmmcopy_prepare
if [ ! -f "genome.fa.map1.bw" ]; then
    /storage/yuhongtaoLab/liangcai/WGS/hmmcopy_utils/util/mappability/generateMap.pl genome/genome.fa -o genome.fa.map1.bw
fi
/storage/yuhongtaoLab/liangcai/WGS/hmmcopy_utils/bin/mapCounter -w 10000 genome.fa.map1.bw > genome.fa.map.w10000.wig
/storage/yuhongtaoLab/liangcai/WGS/hmmcopy_utils/bin/readCounter -w 10000 ${sample_name}_sorted.bam > readcounts.seg.w10000.wig
/storage/yuhongtaoLab/liangcai/WGS/hmmcopy_utils/bin/gcCounter -w 10000 ./genome/genome.fa > genome.gc.w10000.wig

# 4 hmmcopy_plot
source activate
conda deactivate
conda activate alignment
Rscript codes/HMMcopy.R readcounts.seg.w10000.wig genome.gc.w10000.wig genome.fa.map.w10000.wig
