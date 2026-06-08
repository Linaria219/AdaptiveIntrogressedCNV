#!/bin/bash
#SBATCH --job-name=testGSfor40samples
#SBATCH --partition=YOUR_name_of_partition
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH -o testGSfor40samples.out
#SBATCH -e testGSfor40samples1.err

#module load apps/anaconda3/5.2.0
#conda activate genomeSTRiP

which java > /dev/null || exit 1
which Rscript > /dev/null || exit 1
which samtools > /dev/null || exit 1

export SV_DIR="path/to/svtoolkit"
export PATH=${SV_DIR}/bwa:${PATH}
export LD_LIBRARY_PATH=${SV_DIR}/bwa:${LD_LIBRARY_PATH}

classpath="${SV_DIR}/lib/SVToolkit.jar:${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar:${SV_DIR}/lib/gatk/Queue.jar"
rmd="path/to/refmetadata"
rundir="path/to/test_random_40individuals"

java -Xmx16g \
	-cp path/to/svtoolkit/lib/SVToolkit.jar:path/to/svtoolkit/lib/gatk/GenomeAnalysisTK.jar:path/to/svtoolkit/lib/gatk/Queue.jar \
	org.broadinstitute.gatk.queue.QCommandLine \
	-cp path/to/svtoolkit/lib/SVToolkit.jar:path/to/svtoolkit/lib/gatk/GenomeAnalysisTK.jar:path/to/svtoolkit/lib/gatk/Queue.jar \
	-S ${SV_DIR}/qscript/SVPreprocess.q \
	-S ${SV_DIR}/qscript/SVQScript.q \
	-gatk ${SV_DIR}/lib/gatk/GenomeAnalysisTK.jar \
	-configFile ${SV_DIR}/conf/genstrip_parameters.txt \
	-R ${rmd}/Homo_sapiens_assembly38.fasta \
	-genderMapFile ${rundir}/2183_gender.txt \
	-readDepthMaskFile ${rmd}/Homo_sapiens_assembly38.rdmask.bed \
	-genomeMaskFile ${rmd}/Homo_sapiens_assembly38.mask.101.fasta \
	-copyNumberMaskFile ${rmd}/Homo_sapiens_assembly38.gcmask.fasta \
	-ploidyMapFile ${rmd}/Homo_sapiens_assembly38.ploidymap.txt \
	-genderMaskBedFile ${rmd}/Homo_sapiens_assembly38.gendermask.bed \
	-md ${rundir}/metadata \
	-I ${rundir}/20random.list \
	-tempDir ${rundir}/temp \
	-runDirectory ${rundir} \
	-jobLogDir ${rundir}/svpreprocess_joblog \
	-bamFilesAreDisjoint true \
	-computeReadCounts true \
	-P align.altAwareAlignments:true \
	-run

#	-parallelJobs 20 \
#	-run

#	-L chr1 \
#	-L chr2 \
#	-L chr3 \
#	-L chr4 \
#	-L chr5 \
#	-L chr6 \
#	-L chr7 \
#	-L chr8 \
#	-L chr9 \
#	-L chr10 \
#	-L chr11 \
#	-L chr12 \
#	-L chr13 \
#	-L chr14 \
#	-L chr15 \
#	-L chr16 \
#	-L chr17 \
#	-L chr18 \
#	-L chr19 \
#	-L chr20 \
#	-L chr21 \
#	-L chr22 \
#	-L chrX \
#	-L chrY \


