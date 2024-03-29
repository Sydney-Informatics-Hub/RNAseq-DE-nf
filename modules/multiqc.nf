// Define the process
process multiqc {
	// Define directives 
	// See: https://www.nextflow.io/docs/latest/process.html#directives
	debug = true //turn to false to stop printing command stdout to screen
	publishDir "${params.outDir}/multiqc_raw", mode: 'copy'

	// Define input 
	// See: https://www.nextflow.io/docs/latest/process.html#inputs
	input:
	path ('*')

	// Define code to execute 
	// See: https://www.nextflow.io/docs/latest/process.html#script
	script:
	"""
	multiqc .

	"""
 }
