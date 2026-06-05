#!/bin/bash

input_file="/work/home/caizonglin/delly0.7.2/sampleANDpath.txt"
output_dir="/work/home/caizonglin/delly0.7.2/sample1_10"
ref="/work/home/caizonglin/backup_info/GRCh38_full_analysis_set_plus_decoy_hla.fa"

while IFS=$'\t' read -r sample bam_path; do
  cat << EOF > ${output_dir}/${sample}_delly.slurm
#!/bin/bash
#SBATCH --job-name=${sample}_delly
#SBATCH --partition=wzhcnormal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH -o ${output_dir}/${sample}_delly_n1c1.out
#SBATCH -e ${output_dir}/${sample}_delly_n1c1.err

ref="/work/home/caizonglin/backup_info/GRCh38_full_analysis_set_plus_decoy_hla.fa"

delly -t DEL -o ${sample}.del.vcf -g ${ref} -x human.hg38.excl.tsv ${bam_path}
delly -t DUP -o ${sample}.dup.vcf -g ${ref} -x human.hg38.excl.tsv ${bam_path}

bgzip ${sample}.del.vcf
bgzip ${sample}.dup.vcf

tabix -p vcf ${sample}.del.vcf.gz
tabix -p vcf ${sample}.dup.vcf.gz

bcftools concat -a ${sample}.del.vcf.gz ${sample}.dup.vcf.gz -o ${sample}.vcf.gz
EOF
done < ${input_file}

