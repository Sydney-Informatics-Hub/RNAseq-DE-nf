process runTximportCount {
    
    debug = true //turn to false to stop printing command stdout to screen
    publishDir "${params.outDir}/${sampleID}/salmon/converted_counts", mode: 'copy'   
    
    input:
    tuple val(sampleID), path(salmon)
    path makeTx2geneFile

    output:
    path "*gene_tpm.tsv"                 , emit: tpm_gene
    path "*gene_counts.tsv"              , emit: counts_gene
    path "*gene_counts_length_scaled.tsv", emit: counts_gene_length_scaled
    path "*gene_counts_scaled.tsv"       , emit: counts_gene_scaled
    path "*gene_lengths.tsv"             , emit: lengths_gene
    path "*transcript_tpm.tsv"           , emit: tpm_transcript
    path "*transcript_counts.tsv"        , emit: counts_transcript
    path "*transcript_lengths.tsv"       , emit: lengths_transcript
    

    script: // This script is bundled with the pipeline in /bin. 
    
    """
    tximport.r \
        NULL \
        ${salmon} \
        ${sampleID} \
        salmon \
        ${makeTx2geneFile}
    
    """
}