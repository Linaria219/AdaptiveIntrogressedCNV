#!/bin/bash
#SBATCH --job-name=delly1_10
#SBATCH --partition=YOUR_name_of_partition
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH -o aaadelly1_10_n1c1.out
#SBATCH -e aaadelly1_10_n1c1.err
#SBATCH -a 1-10

#conda activate delly0.7.2
#module load apps/samtools/1.9/gcc-7.3.1
#module load bcftools-1.10.2-none
#arrayNUM: 1-10, 11-535, 536-1060, 1061-1585, 1586-2110

ref="path/to/reference.fa"
pathfile="path/to/sampleANDpath.txt"

sample_line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${pathfile})
sample=$(echo $sample_line | awk '{print $1}')
bam_dir=$(echo $sample_line | awk '{print $2}')
echo "${sample}"
echo "${bam_dir}"

delly -t DEL -o ${sample}.del.vcf -g ${ref} -x human.hg38.excl.tsv ${bam_dir}
delly -t DUP -o ${sample}.dup.vcf -g ${ref} -x human.hg38.excl.tsv ${bam_dir}

bgzip ${sample}.del.vcf
bgzip ${sample}.dup.vcf

tabix -p vcf ${sample}.del.vcf.gz
tabix -p vcf ${sample}.dup.vcf.gz

# merge using bcftools or mergeSVcallers
bcftools concat -a ${sample}.del.vcf.gz ${sample}.dup.vcf.gz -o ${sample}.vcf.gz
