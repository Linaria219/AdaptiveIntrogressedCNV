#!/bin/bash

input_file="filename"
output_dir="/path/to/output"
ref="/path/to/reference.fa"

while IFS=$'\t' read -r sample bam_path; do
  cat << EOF > ${output_dir}/${sample}_lumpy.slurm
#!/bin/bash
#SBATCH --job-name=${sample}_lumpy
#SBATCH --partition=YOUR_name_of_partition
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH -o ${sample}_lumpy_n1c20.out
#SBATCH -e ${sample}_lumpy_n1c20.err

samtools sort -@ 15 -m 2G -n -o ${sample}.sort_n.sam ${bam_path}
samblaster -i ${sample}.sort_n.sam --ignoreUnmated --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 | samtools view -hb -o ${sample}.samblaster.bam

samtools view -@ 15 -b -F 1294 ${sample}.samblaster.bam > ${sample}.discordants.unsorted.bam
samtools view -h ${sample}.samblaster.bam | python Python/extractSplitReads_BwaMem -i stdin | samtools view -@ 15 -Sb - > ${sample}.splitters.unsorted.bam
python Python/extractSplitReads_BwaMem -i ${sample}.blaster.bam | samtools view -@ 15 -Sb - > ${sample}.splitters.unsorted.bam

samtools sort -@ 15 -m 2G -o ${sample}.discordants.bam ${sample}.discordants.unsorted.bam
samtools sort -@ 15 -m 2G -o ${sample}.splitters.bam ${sample}.splitters.unsorted.bam

lumpyexpress -x Python/exclude4lumpy.bed -B ${bam_path} -S ${sample}.splitters.bam -D ${sample}.discordants.bam -o ${sample}.vcf
bgzip ${sample}.vcf

rm ${sample}.sort_n.sam
rm ${sample}.samblaster.bam
rm ${sample}.discordants.unsorted.bam
rm ${sample}.splitters.unsorted.bam

EOF
done < ${input_file}

