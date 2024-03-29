#!/bin/bash
  
#PBS -P er01 
#PBS -N RNA-Seq-nf
#PBS -l walltime=01:00:00
#PBS -l ncpus=1
#PBS -l mem=4GB
#PBS -W umask=022
#PBS -q copyq
#PBS -e RNA-Seq-nf.e
#PBS -o RNA-Seq-nf.o
#PBS -l wd
#PBS -l storage=scratch/er01+gdata/er01
#PBS -l jobfs=4GB

# Load singularity and nextflow modules
module purge
#module load nextflow/22.04.3
module load nextflow/23.04.1  
module load singularity
#module load star/2.7.3a

export NXF_SINGULARITY_CACHEDIR=/scratch/$PROJECT/$(whoami)/singularity

# Fill in these variables for your run
#samples=/scratch/er01/ndes8648/pipeline_work/nextflow/INFRA-121-RNASeq-DE/RNAseq-DE-nf/sampleSheet_both_single_paired.csv
#samples=/scratch/er01/ndes8648/pipeline_work/nextflow/INFRA-121-RNASeq-DE/RNAseq-DE-nf/sampleSheet_pairedWithLanes.csv
#samples=/scratch/er01/ndes8648/pipeline_work/nextflow/INFRA-121-RNASeq-DE/RNAseq-DE-nf/sampleSheet_single.csv
samples=/scratch/er01/ndes8648/pipeline_work/nextflow/INFRA-121-RNASeq-DE/RNAseq-DE-nf/sampleSheet_all_paired.csv


#ref=/g/data/er01/SIH-HPC-WGS/Reference/hs38DH.fasta
#dict=/g/data/er01/SIH-HPC-WGS/Reference/hs38DH.dict
#STARRefIndexPath=/g/data/er01/SIH-Gadi-RNAseq/Reference/GRCh38
refFasta=/scratch/er01/ndes8648/pipeline_work/nextflow/INFRA-121-RNASeq-DE/test_data_2024/chromosome22.fasta
#STARRefIndexPath=/scratch/er01/ndes8648/pipeline_work/nextflow/INFRA-121-RNASeq-DE/test_data_2024/STAR_index_singularity_star_2.7.11a
refGtf=/scratch/er01/ndes8648/pipeline_work/nextflow/INFRA-121-RNASeq-DE/test_data_2024/chromosome22.gtf

#SalmonRefIndexPath=
outDir=results_all_paired_UHR_Rep1L1L2_for_PCA
adapters_bbmap=/scratch/er01/ndes8648/pipeline_work/nextflow/INFRA-121-RNASeq-DE/test_data_2024/adapters.fa
transcriptFasta=/scratch/er01/cl9310/1_project/transcriptom/transcriptome.fa
libType=A
#readlen=
strand="reverse"
samples_info=/scratch/er01/ndes8648/pipeline_work/nextflow/INFRA-121-RNASeq-DE/iterations/090224_PCA/tutorial/data-carpentry-rnaseq/data/sample_info_actual_data_all_samples.tsv


# Run the pipeline 
nextflow run main.nf -resume \
        --input ${samples} \
        -profile gadi \
        --whoami $(whoami) --gadi_account $PROJECT \
        --adapters_bbmap ${adapters_bbmap} \
	--refFasta ${refFasta} \
	--refGtf ${refGtf} \
	--transcriptFasta ${transcriptFasta} \
        --libType ${libType} \
        --strand ${strand} \
	--outDir ${outDir} \
        --samples_info ${samples_info}
        
