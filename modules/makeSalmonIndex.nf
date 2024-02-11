// Define the process

process makeSalmonIndex {
 
 debug = true //turn to false to stop printing command stdout to screen
    publishDir "${params.outDir}/INDEX/salmon", mode: 'copy'   
	
    input:
    path refFasta
    path transcriptFasta
    val(NCPUS)

    output:
    path "salmonIndex"
    
    when:
    !file "${params.outDir}/INDEX/salmon".exists()

    script:
    
    """
    grep "^>" <(gunzip -c ${refFasta}) | cut -d " " -f 1 > decoys.txt
    sed -i.bak -e 's/>//g' decoys.txt
    cat ${transcriptFasta} ${refFasta} > gentrome.fasta

    salmon index \\
        --threads ${NCPUS} \\
        -t gentrome.fasta \\
        -d decoys.txt \\
        -i salmon_index \\
        --gencode

    """
}