process makeTx2geneFile {
    
    debug = true //turn to false to stop printing command stdout to screen
    publishDir "${params.outDir}/INDEX", mode: 'copy'

    input:
    path refGtf

    output:
    path 'transcript_to_gene_mapping.txt'

    
    script: // To generate transcript_to_gene_mapping.txt for counts matrix
    """
    gtfToGenePred -genePredExt ${params.refGtf} stdout | cut -f 1,12 > transcript_to_gene_mapping.txt
    
    """
}
