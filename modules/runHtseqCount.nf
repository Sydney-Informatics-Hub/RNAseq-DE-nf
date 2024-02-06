process runHtseqCount {

        // where to publish the outputs
        tag "$uniqueSampleID runHtseqCount"
        publishDir "${params.outDir}/STAR/${uniqueSampleID}", mode:'copy'

        input:
			val uniqueSampleID			
			path uniqueSampleID_bam
			path refGtf
			val(strand)



	output:
		path ("${uniqueSampleID}.counts")

        script:
	
		"""
		htseq-count -f bam -r pos --mode=union -s ${strand} ${uniqueSampleID_bam} ${refGtf} > ${uniqueSampleID}.counts	

		"""
	
	}

	
