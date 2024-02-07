// Define the process
process fastqc {	

    debug = true //turn to false to stop printing command stdout to screen
    publishDir "${params.outDir}/${sampleID}/FastQC", mode: 'copy'

    input: 
    tuple val(sampleID), val(Lane), path(R1), path(R2), val(SEQUENCING_CENTR), val(PLATFORM), val(RUN_TYPE_SINGLE_PAIRED), val(LIBRARY)
    

    output:
    path("${sampleID}*fastqc.html")
    path("${sampleID}*fastqc.zip")

    script:
    """
    mkdir ${sampleID}

    fastqc ${R1} ${R2}

    """
}