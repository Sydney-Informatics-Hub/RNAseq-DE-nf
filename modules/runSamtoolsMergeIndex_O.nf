process runSamtoolsMergeIndex {
  
        // where to publish the outputs
        tag "$sampleID runSTARAlign"
        publishDir "${params.outDir}/STAR/${sampleID}_${lane}", mode:'copy'

        input:
                tuple val(sampleID) , file(R1_fastq), file(R2_fastq), val(seqcentre), val(platform), val(RUN_TYPE_SINGLE_PAIRED), val(lane),val(library)
                        path sampleID_bam
                        val(NCPUS)

        output:
                path ("${sampleID}_${lane}_Aligned.sortedByCoord.out.bam.bai")




        shell:
        '''


        # Merge options to be picked from https://github.com/Sydney-Informatics-Hub/RNASeq-DE/blob/master/Scripts/samtools_merge_index.sh
        # Index for now


        samtools index -@ !{NCPUS} !{sampleID_bam} -o !{sampleID}_!{lane}_Aligned.sortedByCoord.out.bam.bai




        '''
}

