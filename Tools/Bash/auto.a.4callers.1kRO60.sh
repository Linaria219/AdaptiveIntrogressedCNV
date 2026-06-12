#!/bin/bash
ref="path/to/reference/genome.fa"
df="path/to/delly.vcf.gz"
gf="path/to/genomestrip.vcf.gz"
lf="path/to/lumpy.vcf.gz"
wf="path/to/whamg.vcf.gz"
outdir="path/to/output"
sup="path/to/supplement.vcf"

path/to/mergeSVcallers -a ${ref} \
	-f ${df},${gf},${lf},${wf} \
	-t delly,gs,lumpy,whamg \
	-s 1000 \
	-r 0.6 2> ${outdir}/merge4caller.coed.err > ${outdir}/merge4caller.coed.vcf

awk -v insert="$sup" '
BEGIN {
    inserted = 0
}
/^#CHROM/ && !inserted {
    while ((getline line < insert) > 0)
        print line
    close(insert)
    inserted = 1
}
{ print }
' ${outdir}/merge4caller.coed.vcf > ${outdir}/merge4caller.coed.sup.vcf
rm ${outdir}/merge4caller.coed.vcf
bcftools sort -o ${outdir}/merge4caller.coed.srt.vcf ${outdir}/merge4caller.coed.sup.vcf
rm ${outdir}/merge4caller.coed.sup.vcf
bgzip ${outdir}/merge4caller.coed.srt.vcf
tabix -p vcf ${outdir}/merge4caller.coed.srt.vcf.gz

######################### One iteration of merging was done to avoid collapsing unique alleles into the same call.########################################
path/to/mergeSVcallers -a ${ref} \
	-f ${outdir}/merge4caller.coed.srt.vcf.gz \
	-t iteration1 \
	-s 1000 \
	-r 0.6 2> ${outdir}/merge4caller.coed.iteration1.err > ${outdir}/merge4caller.coed.iteration1.vcf

awk -v insert="$sup" '
BEGIN {
    inserted = 0
}
/^#CHROM/ && !inserted {
    while ((getline line < insert) > 0)
        print line
    close(insert)
    inserted = 1
}
{ print }
' ${outdir}/merge4caller.coed.iteration1.vcf > ${outdir}/merge4caller.coed.iteration1.sup.vcf
rm ${outdir}/merge4caller.coed.iteration1.vcf
bcftools sort -o ${outdir}/merge4caller.coed.iteration1.srt.vcf ${outdir}/merge4caller.coed.iteration1.sup.vcf
rm ${outdir}/merge4caller.coed.iteration1.sup.vcf
bgzip ${outdir}/merge4caller.coed.iteration1.srt.vcf
tabix -p vcf ${outdir}/merge4caller.coed.iteration1.srt.vcf.gz

