process runSTARAlign {

        // where to publish the outputs
        tag "$sampleID runSTARAlign"
        publishDir "${params.outDir}/${sampleID}/STAR", mode:'copy'

        input:
        	tuple val(sampleID) , val(lane), file(R1_fastq), file(R2_fastq), val(seqcentre), val(platform), val(RUN_TYPE_SINGLE_PAIRED) ,val(library)
			val(NCPUS)	
			path STARRefIndexPath
			path R1_trimmed
			path R2_trimmed

        output:
           path ("${sampleID}_${lane}_Aligned.sortedByCoord.out.bam")
	   path ("${sampleID}_${lane}_SJ.out.tab")

        script:
	
		"""

		flowcell=1
                #lane=1
	
		if [ "${RUN_TYPE_SINGLE_PAIRED}" == 'PAIRED' ]; then


		# Mapping
		STAR \
            		--runThreadN ${NCPUS} \
            		--outBAMsortingThreadN ${NCPUS} \
            		--genomeDir ${STARRefIndexPath} \
            		--outBAMsortingBinsN 100 \
            		--quantMode GeneCounts \
            		--readFilesCommand zcat \
            		--readFilesIn ${R1_trimmed} ${R2_trimmed} \
            		--outSAMattrRGline ID:flowcell.${lane} PU:flowcell.${lane}.${sampleID} SM:${sampleID} PL:${platform} CN:${seqcentre} LB:${library} \
            		--outSAMtype BAM SortedByCoordinate \
            		--outReadsUnmapped Fastx \
            		--outSAMunmapped Within KeepPairs \
            		--outFileNamePrefix ${sampleID}_${lane}_

		else

		STAR \
			--runThreadN ${NCPUS} \
			--outBAMsortingThreadN ${NCPUS} \
			--genomeDir ${STARRefIndexPath} \
			--quantMode GeneCounts \
        		--outBAMsortingBinsN 100 \
			--readFilesCommand zcat \
			--readFilesIn ${R1_trimmed} \
			--outSAMattrRGline ID:flowcell.${lane} PU:flowcell.${lane}.${sampleID} SM:${sampleID} PL:${platform} CN:${seqcentre} LB:${library} \
			--outSAMtype BAM SortedByCoordinate \
			--outReadsUnmapped Fastx \
			--outFileNamePrefix ${sampleID}_${lane}_
		

		fi
			



		"""
	
	}

	
