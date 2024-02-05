// Define the process
process checkCohort {	
	
// check the existence of input files 
input:
	path input

output:
	path "samples.txt"
		
	script:
	"""
	cat "${params.input}" > samples.txt
	"""
}