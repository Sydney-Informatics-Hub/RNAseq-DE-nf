#!/usr/bin/env nextflow

// =================================================================
// main.nf is the pipeline script for a nextflow pipeline
// Should contain the following sections:
	// Process definitions
    // Channel definitions
    // Workflow structure
	// Workflow summary logs 

// Examples are included for each section. Remove them and replace
// with project-specific code. For more information see:
// https://www.nextflow.io/docs/latest/index.html.
//
// ====================================================================

// Import processes or subworkflows to be run in the workflow
// Each of these is a separate .nf script saved in modules/ directory
// See https://training.nextflow.io/basic_training/modules/#importing-modules 
include { checkCohort} from './modules/checkCohort.nf'
include { fastqc     } from './modules/fastqc.nf'
include { multiqc    } from './modules/multiqc.nf'
include { bbduk      } from './modules/bbduk.nf'
include { makeSTARIndex                                  } from './modules/makeSTARIndex'
include { runSTARAlign                                   } from './modules/runSTARAlign'
include { runSamtoolsMergeIndex                          } from './modules/runSamtoolsMergeIndex'
include { runHtseqCount                                   } from './modules/runHtseqCount'
include { mergeHtseqCounts                                } from './modules/mergeHtseqCounts'
include { makeSalmonIndex                                } from './modules/makeSalmonIndex.nf'
include { runSalmonAlign                                } from './modules/runSalmonAlign.nf'
include { createPCAFromCounts                            } from './modules/createPCAFromCounts'

include { convertGtfToBED                                } from './modules/convertGtfToBED'
include { getMappingMetricRSeQC                          } from './modules/getMappingMetricRSeQC'

// Print a header for your pipeline 
log.info """\

====================================================================
RNASEQ DE NF
====================================================================


Created by the Sydney Informatics Hub, University of Sydney
Find documentation @ https://github.com/Sydney-Informatics-Hub/RNAseq-DE-nf
Cite this pipeline @ INSERT DOI

====================================================================
Workflow run parameters 
====================================================================
input       : ${params.input}
outDir      : ${params.outDir}
workDir     : ${workflow.workDir}
====================================================================

"""

/// Help function 
// This is an example of how to set out the help function that 
// will be run if run command is incorrect or missing. 

def helpMessage() {
    log.info"""
  Usage:   nextflow run main.nf --input samples.csv \
                     			--refFasta /path/to/ref.fasta --refGtf /path/to/ref.gtf

  Required Arguments:
    --input                                   Full path and name of sample input file (csv format)
    --refFasta                                Full path and name of reference genome (fasta format)
    --refGtf								  Full path and name of reference annotation (gtf format)

  Optional Arguments:
    --outDir                        Specify name of results directory.


  HPC accounting arguments:
    --whoami                    HPC user name (Setonix or Gadi HPC)
    --gadi_account              Project accounting code for NCI Gadi (e.g. aa00)
  """.stripIndent()
}


// Define workflow structure. Include some input/runtime tests here.
// See https://www.nextflow.io/docs/latest/dsl2.html?highlight=workflow#workflow
workflow {

// Show help message if --help is run or (||) a required parameter (input) is not provided

if ( params.help || params.input == false ){   
// Invoke the help function above and exit
	helpMessage()
	exit 1
	// consider adding some extra contigencies here.
	// could validate path of all input files in list?
	// could validate indexes for reference exist?

// If none of the above are a problem, then run the workflow
} else {
	
// Define channels 
// See https://www.nextflow.io/docs/latest/channel.html#channels
// See https://training.nextflow.io/basic_training/channels/
// check the existence of input files  
	checkCohort(Channel.fromPath(params.input, checkIfExists: true))

// make a input channel and map all tuples

inputs = checkCohort.out
                .splitCsv(header: true, sep:",")
                .map { row ->
                def R2File = row.R2 ? file(row.R2) : file(emptyFilePath) // Provide a default empty file path
                tuple(row.sampleID, row.Lane, file(row.R1), R2File, row.SEQUENCING_CENTRE, row.PLATFORM, row.RUN_TYPE_SINGLE_PAIRED, row.LIBRARY)
                        }

// Run fastqc
// See https://training.nextflow.io/basic_training/processes/#inputs 
fastqc(inputs)

	//multiqc(fastqc.out.collect())

// Run trimming 
bbduk(params.adapters_bbmap, inputs)

// Define star/salmon input
align_input = bbduk.out.trimmed_fq
  //.view() //ADDED THIS TO VISUALISE OUT STRUCTURE FOR DEBUGGING PURPOSES
  .map { tuple ->
  // Extracting values and paths from the tuple produced by bbduk.out.trimmed_fq
  def sampleID = tuple[0]
  def lane = tuple[1]
  def runType = tuple[2]
  def platform = tuple[3]  
  def sequencingCentre = tuple[4]
  def library = tuple[5]
  def r1Path = tuple[6]
  def r2Path = tuple[7]

    if (runType == "PAIRED") {
        return [sampleID, lane, runType, platform, sequencingCentre, library, r1Path, r2Path]
    } else if (runType == "SINGLE") {
        // For SINGLE runs, emit a tuple with only R1 path and an empty string for R2
        return [sampleID, lane, runType, platform, sequencingCentre, library, r1Path, ""]
    }
}

alignmentInputSalmon = align_input.groupTuple()

// Run STAR index and alignment

//STAR_ref_index_path = "$PWD/${params.outDir}/INDEX/STAR/STARGeneratedIndexPath"

//if (!file(STAR_ref_index_path).exists()) {
        // Make STAR index and then align
        makeSTARIndex(params.refFasta,params.refGtf)
        runSTARAlign(makeSTARIndex.out.STAR_ref_index_path,align_input)

//} else if (file(STAR_ref_index_path).exists()){
 //       runSTARAlign(STAR_ref_index_path,align_input)
 //       }



samtoolsMergeInput = runSTARAlign.out.sample_lane_bam
  //.view() //ADDED THIS TO VISUALISE OUT STRUCTURE FOR DEBUGGING PURPOSES
  .map { tuple ->
  // Extracting values and paths from the tuple produced by runSTARAlign.out.sample_lane_bam
  def sampleID = tuple[0]
  def path_bam = tuple[1]
  def path_SJ = tuple[2]


  return [sampleID, path_bam, path_SJ]

} | groupTuple()


// Merge lane-bams and Index final bam
runSamtoolsMergeIndex(samtoolsMergeInput)


merged_input=runSamtoolsMergeIndex.out.final_bam
  //.view()
  .map { tuple ->
  // Extracting values and paths from the tuple produced by runSamtoolsMergeIndex.out.final_bam
  def sampleID = tuple[0]
  def pathBam = tuple[1]
  def pathBamBai = tuple[2]

  return [sampleID, pathBam, pathBamBai]
}

// Run HTSeq-Count
runHtseqCount(merged_input,params.refGtf,params.strand)

mergeHtseqCounts(runHtseqCount.out.sampleIDCounts.collect())


// Run Salmon Index and alignment

        makeSalmonIndex(params.refFasta,params.transcriptFasta)
        runSalmonAlign(makeSalmonIndex.out,params.libType,align_input)
        
// Create PCA 
createPCAFromCounts(mergeHtseqCounts.out.merged_counts_STAR,params.samples_info)

// Create BED and use ot for STAR metrics
convertGtfToBED(params.refGtf)
getMappingMetricRSeQC(merged_input,convertGtfToBED.out)

}}

// Print workflow execution summary 
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