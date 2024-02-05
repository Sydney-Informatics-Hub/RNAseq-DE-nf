// Define the process
process fastqc {	

    debug = true //turn to false to stop printing command stdout to screen
    publishDir "${params.outDir}/${sampleID}", mode: 'copy'

    input: 
    tuple val(sampleID), path(R1), path(R2)

    output:
    path("${sampleID}*fastqc.html")
    path("${sampleID}*fastqc.zip")

    script:
    """
    mkdir ${sampleID}
    fastqc ${R1}${R2}

    """
}