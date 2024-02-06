// Define the process
process multiqc {
	// Define directives 
	// See: https://www.nextflow.io/docs/latest/process.html#directives
	debug = true //turn to false to stop printing command stdout to screen
	publishDir "${params.outDir}", mode: 'copy'

	// Define input 
	// See: https://www.nextflow.io/docs/latest/process.html#inputs
	input:
	path ('*')
	tuple val(sampleID), val(Lane), path(R1), path(R2), val(SEQUENCING_CENTR), val(PLATFORM), val(RUN_TYPE_SINGLE_PAIRED), val(LIBRARY)

	// Define output(s)
	// See: https://www.nextflow.io/docs/latest/process.html#outputs
	output:
	path ('multiqc_report.html')
	path ('multiqc_data')

	// Define code to execute 
	// See: https://www.nextflow.io/docs/latest/process.html#script
	script:
	"""
	multiqc .

	"""
 }