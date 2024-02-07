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


inputs = checkCohort.out
		.splitCsv(header: true, sep:",")
		.map { row -> tuple(row.sampleID, row.Lane, file(row.R1), file(row.R2), row.SEQUENCING_CENTRE, row.PLATFORM, row.RUN_TYPE_SINGLE_PAIRED, row.LIBRARY)}

// Run fastqc
// See https://training.nextflow.io/basic_training/processes/#inputs 
	fastqc(inputs)
	multiqc(fastqc.out[1].collect())
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