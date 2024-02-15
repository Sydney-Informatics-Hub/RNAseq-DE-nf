// Define the process

process runSalmonAlign {
 
 debug = true //turn to false to stop printing command stdout to screen
    publishDir "${params.outDir}/${sampleID}/salmon/quant_${Lane}", mode: 'copy'   
	
    input:
		path salmonIndex
        val(libType)
		tuple val(sampleID), val(Lane) , val(RUN_TYPE_SINGLE_PAIRED), val(platform), val(seqcentre), val(library), path(r1Path), path(r2Path)

    output:
    path salmon
    
    script:
    """

    salmon quant \
        --threads ${task.cpus} \
        -i ${salmonIndex} \
        -l ${libType} \
        -1 ${r1Path} \
        -2 ${r2Path} \
        --validateMappings \
        -o salmon
    """
}