#!/usr/bin/env nextflow

// Import subworkflows to be run in the workflow
include { checkInputs                                    } from './modules/checkCohort'
include { convertGtfToBED                                } from './modules/convertGtfToBED'

include { makeSTARIndex					 } from './modules/makeSTARIndex'

// Mapping is done at lane level (per sample pair as is in sampleSheet) and then merge bams 
include { runSTARAlign                                   } from './modules/runSTARAlign' 
include { runSamtoolsMergeIndex		                 } from './modules/runSamtoolsMergeIndex'

include { runHtseqCount                                   } from './modules/runHtseqCount'
include { mergeHtseqCounts                                } from './modules/mergeHtseqCounts'
include { getMappingMetricRSeQC                           } from './modules/getMappingMetricRSeQC'


/// Print a header for your pipeline 

log.info """\

===================================================================
===================================================================
SOMATIC SHORT V - NF 
===================================================================
===================================================================

Created by the Sydney Informatics Hub, University of Sydney

Documentation	@ https://github.com/Sydney-Informatics-Hub/Somatic-shortV-nf

Cite					@ https://doi.org/10.48546/workflowhub.workflow.691.1

Log issues    @ https://github.com/Sydney-Informatics-Hub/Somatic-shortV-nf/issues

All the default parameters are set in `nextflow.config`

=======================================================================================
Workflow run parameters 
=======================================================================================
version                    : ${params.version}
input                      : ${params.input}
reference                  : ${params.ref}
dict                       : ${params.dict} 
common_biallelic_variants  : ${params.common_biallelic_variants}
ponvcf                     : ${params.ponvcf}
outDir                     : ${params.outDir}
workDir                    : ${workflow.workDir}

=======================================================================================

 """

/// Help function 
// This is an example of how to set out the help function that 
// will be run if run command is incorrect (if set in workflow) 
// or missing/  

def helpMessage() {
    log.info"""
  Usage:   nextflow run main.nf --input samples.csv \
                     --ref /path/to/ref.fasta --dict /path/to/ref.dict \
                     --ponvcf /path/to/pon \
                     --common_biallelic_variants /path/to/common_biallelic_variants

  Required Arguments:
    --input		                      Full path and name of sample input file (csv format)
    --ref			                      Full path and name of reference genome (fasta format)
    --dict                          Full path and name of reference genome dictionary file (dict format)
    --ponvcf                        Full path and name of the Panel of Normals file (vcf format)
    --common_biallelic_variants     Full path and name of the common biallelic variant resources file (vcf format)
	
  Optional Arguments:
    --outDir                        Specify name of results directory. 
    --number_of_intervals           Define a specific number genomic-intervals for parallelisation


  HPC accounting arguments:
    --whoami                    HPC user name (Setonix or Gadi HPC)
    --gadi_account              Project accounting code for NCI Gadi (e.g. aa00)
  """.stripIndent()
}

/// Main workflow structure. 

workflow {

// Show help message if --help is run or if any required params are not 
// provided at runtime

  if ( params.help == true || params.input == false)
	{   

          // Invoke the help function above and exit
          helpMessage()
          exit 1
        
	} 


	else 
	{
	
  // Check inputs file exists

	//Input samplesheet looks like this
	//ID1,ID1_L1_R1.fastq,ID1_L1_R2.fastq
	//ID1,ID1_L2_R1.fastq,ID1_L2_R2.fastq
  //ID2,ID2_L1_R1.fastq,ID2_L1_R2.fastq
	//ID2,ID2_L2_R1.fastq,ID2_L2_R2.fastq

	checkInputs(Channel.fromPath(params.input, checkIfExists: true))


	uniqueSampleIDs = checkInputs.out.flatMap { filePath ->
		    file(filePath).text.readLines().drop(1).collect { line -> line.split(',')[0] }
						}
			.distinct()
	
	uniqueSampleIDsList = uniqueSampleIDs.toList()

	uniqueSampleIDsList.view()


	// Split cohort file to collect info for each sample
  //fastq_eachLane_ch = checkInputs.out
  //              .splitCsv(header: true, sep:",")
  //              .map { row -> tuple(row.sampleID, file(row.R1), file(row.R2),row.SEQUENCING_CENTRE,row.PLATFORM,row.RUN_TYPE_SINGLE_PAIRED,row.Lane,row.LIBRARY)}

// Define a valid empty file path using $PWD
def emptyFilePath = "$PWD/empty_file.txt"

// Check if the empty file exists, create it if necessary
if (!file(emptyFilePath).exists()) {
    file(emptyFilePath).text = ""
}


fastq_eachLane_ch = checkInputs.out
	    	.splitCsv(header: true, sep:",")
    		.map { row ->
        	def R2File = row.R2 ? file(row.R2) : file(emptyFilePath) // Provide a default empty file path
        	tuple(row.sampleID, file(row.R1), R2File, row.SEQUENCING_CENTRE, row.PLATFORM, row.RUN_TYPE_SINGLE_PAIRED, row.Lane, row.LIBRARY)
    			}




//ACTUALL - It is required to map per lane level and then merge !
//	=========== >Merge lane level to sample level BAMs: 

//UHR_Rep1,1,/scratch/er01/ndes8648/pipeline_work/nextflow/INFRA-121-RNASeq-DE/temp_git_repo/test_data_2024/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz,/scratch/er01/ndes8648/pipeline_work/nextflow/INFRA-121-RNASeq-DE/temp_git_repo/test_data_2024/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz,KCCG,ILLUMINA,PAIRED,1
//HBR_Rep3,1,/scratch/er01/ndes8648/pipeline_work/nextflow/INFRA-121-RNASeq-DE/temp_git_repo/test_data_2024/HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz,/scratch/er01/ndes8648/pipeline_work/nextflow/INFRA-121-RNASeq-DE/temp_git_repo/test_data_2024/HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz,KCCG,ILLUMINA,PAIRED,1
//HBR_Rep2,1,/scratch/er01/ndes8648/pipeline_work/nextflow/INFRA-121-RNASeq-DE/temp_git_repo/test_data_2024/HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz,/scratch/er01/ndes8648/pipeline_work/nextflow/INFRA-121-RNASeq-DE/temp_git_repo/test_data_2024/HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz,KCCG,ILLUMINA,PAIRED,1
//HBR_Rep1,1,/scratch/er01/ndes8648/pipeline_work/nextflow/INFRA-121-RNASeq-DE/temp_git_repo/test_data_2024/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz,/scratch/er01/ndes8648/pipeline_work/nextflow/INFRA-121-RNASeq-DE/temp_git_repo/test_data_2024/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz,KCCG,ILLUMINA,PAIRED,1


//Run the processes 


convertGtfToBED(params.refGtf)

makeSTARIndex(params.refFasta,params.refGtf,params.NCPUS)


// STAR Alignment

// (A) Use available STAR index
//runSTARAlign(fastq_eachLane_ch,runBbduk_trim.out[0],runBbduk_trim.out[1],params.NCPUS,params.STARRefIndexPath)

// (B) Align
runSTARAlign(fastq_eachLane_ch,params.NCPUS,makeSTARIndex.out)



// Merge and Index

runSamtoolsMergeIndex(uniqueSampleIDs,runSTARAlign.out[0].collect(),params.NCPUS)


//getMappingMetricRSeQC(fastq_eachLane_ch,convertGtfToBED.out,runSTARAlign.out[0],runSamtoolsMergeIndex.out)
getMappingMetricRSeQC(runSamtoolsMergeIndex.out[2],convertGtfToBED.out,runSamtoolsMergeIndex.out[0],runSamtoolsMergeIndex.out[1])


// HTseqCount
runHtseqCount(runSamtoolsMergeIndex.out[2],runSamtoolsMergeIndex.out[0],params.refGtf,params.strand)
mergeHtseqCounts(runHtseqCount.out.collect())


}}

workflow.onComplete {
  summary = """
=======================================================================================
Workflow execution summary
=======================================================================================

Duration    : ${workflow.duration}
Success     : ${workflow.success}
workDir     : ${workflow.workDir}
Exit status : ${workflow.exitStatus}
outDir      : ${params.outDir}

=======================================================================================
  """
  println summary

}

