// Define the process
process fastqc {	

    debug = true //turn to false to stop printing command stdout to screen
    publishDir "${params.output}/${sampleID}", mode: 'copy'

    input: 
    tuple val(sampleID), path(fastq)

    output:
    path("${sampleID}*fastqc.html")
    path("${sampleID}*fastqc.zip")

    script:
    """
    mkdir ${sampleID}
    fastqc ${fastq}

    """
}