process runSamtoolsMergeIndex {

	// where to publish the outputs
        tag "$sampleID runSamtoolsMergeIndex"
        publishDir "${params.outDir}/STAR/${uniqueSampleID}", mode:'copy'

  
        input:
        	each uniqueSampleID        
	        path(sampleID_bams)
		val(NCPUS)

        output:
                path ("${final}.bam.bai")


	script:

	"""

	 # Merge options to be picked from https://github.com/Sydney-Informatics-Hub/RNASeq-DE/blob/master/Scripts/samtools_merge_index.sh

        outDir=${params.outDir}

        sampbams=()
        sampbams+=("$(ls ${outDir}/STAR/${uniqueSampleID}/${uniqueSampleID}*_Aligned.sortedByCoord.out.bam)")
        final=${sampleid}.final.bam

        samtools merge -f -@ ${NCPUS} ${final} ${sampbams[@]}


        # Index for now
        samtools index -@ ${NCPUS} ${final} -o !{final}.bam.bai




	"""

}
