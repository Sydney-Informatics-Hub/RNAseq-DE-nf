// Define the process
process fastqc {	
///	LEAVING DEBUGGING ON WILL PRINT THE STDOUT OF THE COMMAND TO THE SCREEN
    debug = false //turn to false to stop printing command stdout to screen
    tag "$sampleID fastQC"
    publishDir "${params.outDir}/${sampleID}/FastQC", mode: 'copy'

    input: 
    	tuple val(sampleID), val(Lane), path(R1), path(R2), val(SEQUENCING_CENTR), val(PLATFORM), val(RUN_TYPE_SINGLE_PAIRED), val(LIBRARY)

    output:
	path ("*_fastqc.html")

    script:
    """
	if [ "${RUN_TYPE_SINGLE_PAIRED}" == 'PAIRED' ]; then
    		fastqc -t ${task.cpus}  ${R1} ${R2} 
	else
		fastqc -t ${task.cpus}  ${R1}

	fi
	
    """
}
