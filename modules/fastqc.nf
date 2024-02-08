// Define the process
process fastqc {	

    debug = true //turn to false to stop printing command stdout to screen
    tag "$sampleID fastQC"
    publishDir "${params.outDir}/${sampleID}/FastQC", mode: 'copy'

    input: 
    	tuple val(sampleID), val(Lane), path(R1), path(R2), val(SEQUENCING_CENTR), val(PLATFORM), val(RUN_TYPE_SINGLE_PAIRED), val(LIBRARY)
    	val(NCPUS)

    output:
	path ("*_fastqc.html")

    script:
    """
    	mkdir ${sampleID}


	if [ "${RUN_TYPE_SINGLE_PAIRED}" == 'PAIRED' ]; then
    		fastqc -t ${NCPUS}  ${R1} ${R2} 
	else
		fastqc -t ${NCPUS}  ${R1}

	fi
	
    """
}
