// Define the process

process makeSalmonIndex {
 
 debug = true //turn to false to stop printing command stdout to screen
    publishDir "${params.outDir}/INDEX", mode: 'copy'   
	
    input:
    path refFasta
    path transcriptFasta

    output:
    path "salmonIndex"
    
    script:
    
    """
    grep "^>" ${refFasta} | cut -d " " -f 1 > decoys.txt
    sed -i.bak -e 's/>//g' decoys.txt
    cat ${transcriptFasta} ${refFasta} > gentrome.fasta

    salmon index \
        --threads ${task.cpus} \
        -t gentrome.fasta \
        -d decoys.txt \
        -i salmonIndex \
        --gencode

    """
}