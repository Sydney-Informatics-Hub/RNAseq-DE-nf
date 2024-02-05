process runHtseqCount {

        // where to publish the outputs
        tag "$sampleID runSTARAlign"
        publishDir "${params.outDir}/STAR/${sampleID}", mode:'copy'

        input:
        	tuple val(sampleID) , file(R1_fastq), file(R2_fastq), val(seqcentre), val(platform), val(RUN_TYPE_SINGLE_PAIRED), val(lane),val(library)
			path sampleID_bam
			path refGtf
			val(strand)


	output:
		path ("${sampleID}.counts")

        script:
	
		"""
		htseq-count -f bam -r pos --mode=union -s ${strand} ${sampleID_bam} ${refGtf} > ${sampleID}.counts	

		"""
	
	}

	
