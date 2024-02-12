// Define the process

process runSalmonAlign {
 
 debug = true //turn to false to stop printing command stdout to screen
    publishDir "${params.outDir}/${sampleID}/salmon/${Lane}", mode: 'copy'   
	
    input:
        val(NCPUS)	
		path salmonIndex
        val(libType)
		tuple val(sampleID), val(Lane) , val(RUN_TYPE_SINGLE_PAIRED), val(platform), val(seqcentre), val(library), path(sampleID_lane_Trimmed_R1_fastq)
        tuple val(sampleID), val(Lane) , val(RUN_TYPE_SINGLE_PAIRED), val(platform), val(seqcentre), val(library), path(sampleID_lane_Trimmed_R2_fastq)


    output:
    path salmon
    
    script:
    
    """
    salmon quant \
        --threads ${NCPUS} \
        -i ${salmonIndex} \
        -l ${libType} \
        -1 ${sampleID_lane_Trimmed_R1_fastq} \
        -2 ${sampleID_lane_Trimmed_R2_fastq} \
        --validateMappings \
        -o salmon

    """
}

