process runHtseqCount {

        // where to publish the outputs
        tag "$sampleID runHtseqCount"
        publishDir "${params.outDir}/${sampleID}/STAR", mode:'copy'

        input:
			tuple val(sampleID), path(final_bam_file), path(final_bam_index)
			path refGtf
			val strand



	output:
		path ("${sampleID}.counts"), emit: sampleIDCounts

        script:
	
		"""
		 htseq-count -f bam -r pos --mode=union -s ${strand} ${final_bam_file} ${refGtf} > ${sampleID}.counts	

		#touch ${sampleID}.counts
		"""
	
	}

	
