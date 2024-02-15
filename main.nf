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
//include { makeSTARIndex                                  } from './modules/makeSTARIndex'
//include { runSTARAlign                                   } from './modules/runSTARAlign'
//include { runSamtoolsMergeIndex                          } from './modules/runSamtoolsMergeIndex'
//include { runHtseqCount                                   } from './modules/runHtseqCount'
//include { mergeHtseqCounts                                } from './modules/mergeHtseqCounts'
include { makeSalmonIndex                                } from './modules/makeSalmonIndex.nf'
include { runSalmonAlign                                } from './modules/runSalmonAlign.nf'

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

// Create a list of unique sampleIDs
/// HAVING ONLY WORKED ON PASSING BBDUK OUTPUT TO STAR AT THE LANE LEVEL I AM NOT SURE WHY WE NEED THIS
/// UNIQUE SAMPLE IDS ARE DEFINED IN YOUR INPUT CHANNEL FROM YOUR SAMPLESHEET 
/// AT WHAT POINT WOULD YOU NEED UNIQUE SAMPLE IDS THAT ARE DISCONNECTED FROM THEIR R1/2 FILES AND METADATA?
/// IF THIS IS FOR MERGING LANE LEVEL FILES, THEN YOU WOULD NEED TO CREATE A CHANNEL FOR THE BAM MERGING STEP
/// THAT GROUPS BAMS BASED ON THEIR SAMPLEID AND THEN MERGES THEM 

// To account for missing R2 file when single-end 
// Define a valid empty file path using $PWD

// Check if the empty file exists, create it if necessary
/// AGAIN, I'M NOT SURE WHY WE NEED THIS
/// SEE WHAT I DID BELOW LINES 166-171 TO DYNAMICALLY HANDLE OPTIONAL R2 FILE
/// CREATING EMPTY FILES IS GOING TO AFFECT YOUR INODE LIMITS ON THE FILESYSTEM

inputs = checkCohort.out
                .splitCsv(header: true, sep:",")
                .map { row ->
                def R2File = row.R2 ? file(row.R2) : file(emptyFilePath) // Provide a default empty file path
                tuple(row.sampleID, row.Lane, file(row.R1), R2File, row.SEQUENCING_CENTRE, row.PLATFORM, row.RUN_TYPE_SINGLE_PAIRED, row.LIBRARY)
                        }

// Run fastqc
// See https://training.nextflow.io/basic_training/processes/#inputs 
fastqc(inputs)
  /// AS DISCUSSED IT DOESN'T MAKE SENSE TO RUN MULTIQC HERE, AS IT IS NOT THE ONLY QC METRICS BEING CREATED
  /// JUST RUN AT THE END SO YOU CAN COLLECT ALL THE FASTQC OUTPUTS AND ALL OTHER QC METRICS GENERATED THROUGHOUT
	//multiqc(fastqc.out.collect())

// Run trimming 
bbduk(params.adapters_bbmap, inputs)
/// SUGGEST RUNNING FASTQC AGAIN HERE, BUT ON BBDUK OUTPUT 
/// WILL NEED TO CREATE A SPECIFIC BBDUK_FASTQC CHANNEL, SEPARATE FROM STAR OR SALMON
/// IF YOU DON'T, YOU'LL CONSUME BBDUK OUTPUT IN THE MULTIQC STEP, AND YOU WON'T BE ABLE TO USE IT FOR STAR OR SALMON
/// SEE: https://www.nextflow.io/docs/latest/channel.html

// Define star/salmon input
/// THIS IS HOW YOU CAN TAKE THE OUTPUT FROM A PROCESS, AND USE IT AS INPUT TO ANOTHER PROCESS
/// KEEP IN MIND THAT NEXTFLOW IS BASED ON GROOVY
/// SO, WHEN YOU'RE WORKING WITH TUPLES IN NEXTFLOW, YOU'RE WORKING WITH GROOVY LISTS
/// WHEN YOU CAN'T FIND A STRAIGHTFORWARD EXAMPLE IN NEXTFLOW DOCS, LOOK FOR GROOVY EXAMPLES

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

/// THIS WAS TRICKY, TO CAPTURE FLEXIBILITY FOR SINGLE VS PAIRED, NEED TO DEFINE DIFFERENT TUPLES FOR EACH
// If runType is PAIRED, emit a tuple with both R1 and R2 paths
    if (runType == "PAIRED") {
        return [sampleID, lane, runType, platform, sequencingCentre, library, r1Path, r2Path]
    } else if (runType == "SINGLE") {
        // For SINGLE runs, emit a tuple with only R1 path and an empty string for R2
        return [sampleID, lane, runType, platform, sequencingCentre, library, r1Path, ""]
    }
}
| groupTuple
| view

// Run STAR index and alignment
/// THIS LOGIC DOESN'T MAKE SENSE. WHY ARE WE RELIANT ON SOMETHING SAVED TO RESULTS? 
/// SHOULD BE PICKING UP THE INDEX FROM THE REF FASTA DIRECTORY SUPPLIED BY THE USER IF IT EXISTS
/// CREATING TEMPORARY WORKAROUNDS LIKE THIS CREATE MORE WORK FOR YOU IN THE LONG RUN
//STAR_ref_index_path = "$PWD/${params.outDir}/INDEX/STAR/STARGeneratedIndexPath"

//if (!file(STAR_ref_index_path).exists()) {
        // Make STAR index and then align
//        makeSTARIndex(params.refFasta,params.refGtf)
//        runSTARAlign(makeSTARIndex.out.STAR_INDEX,align_input)

//} else if (file(STAR_ref_index_path).exists()){
//        runSTARAlign(makeSTARIndex.out.STAR_ref_index_path,align_input)
//        }

// Merge lane-bams and Index final bam
//runSamtoolsMergeIndex(uniqueSampleIDs,runSTARAlign.out.sampleID_lane_bam.collect(),params.NCPUS)

// Run HTSeq-Count
//runHtseqCount(runSamtoolsMergeIndex.out[2],runSamtoolsMergeIndex.out[0],params.refGtf,params.strand)
//mergeHtseqCounts(runHtseqCount.out.collect())

// Run Salmon Index and alignment

        makeSalmonIndex(params.refFasta,params.transcriptFasta)
        runSalmonAlign(makeSalmonIndex.out,params.libType,align_input)
        
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