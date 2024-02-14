process runSTARAlign {

        // where to publish the outputs
        tag "$sampleID runSTARAlign"
        publishDir "${params.outDir}/${sampleID}/STAR", mode:'copy'

        input:
		val(NCPUS)	
		path STARRefIndexPath
		tuple val(sampleID), val(Lane) , val(RUN_TYPE_SINGLE_PAIRED), val(platform), val(seqcentre), val(library), path(sampleID_lane_Trimmed_R1_fastq)
        	tuple val(sampleID), val(Lane) , val(RUN_TYPE_SINGLE_PAIRED), val(platform), val(seqcentre), val(library), path(sampleID_lane_Trimmed_R2_fastq)


        output:
           	path ("${sampleID}_${Lane}_Aligned.sortedByCoord.out.bam") , emit: sampleID_lane_bam
	   	path ("${sampleID}_${Lane}_SJ.out.tab")			   , emit: sampleID_lane_SJ_tab

        script:
	
		"""


		flowcell=1
	
		if [ "${RUN_TYPE_SINGLE_PAIRED}" == 'PAIRED' ]; then


		# Mapping
		STAR \
            		--runThreadN ${NCPUS} \
            		--outBAMsortingThreadN ${NCPUS} \
            		--genomeDir ${STARRefIndexPath} \
            		--outBAMsortingBinsN 100 \
            		--quantMode GeneCounts \
            		--readFilesCommand zcat \
            		--readFilesIn ${sampleID_lane_Trimmed_R1_fastq} ${sampleID_lane_Trimmed_R2_fastq} \
            		--outSAMattrRGline ID:flowcell.${Lane} PU:flowcell.${Lane}.${sampleID} SM:${sampleID} PL:${platform} CN:${seqcentre} LB:${library} \
            		--outSAMtype BAM SortedByCoordinate \
            		--outReadsUnmapped Fastx \
            		--outSAMunmapped Within KeepPairs \
            		--outFileNamePrefix ${sampleID}_${Lane}_

		else

		STAR \
			--runThreadN ${NCPUS} \
			--outBAMsortingThreadN ${NCPUS} \
			--genomeDir ${STARRefIndexPath} \
			--quantMode GeneCounts \
        		--outBAMsortingBinsN 100 \
			--readFilesCommand zcat \
			--readFilesIn ${sampleID_lane_Trimmed_R1_fastq} \
			--outSAMattrRGline ID:flowcell.${Lane} PU:flowcell.${Lane}.${sampleID} SM:${sampleID} PL:${platform} CN:${seqcentre} LB:${library} \
			--outSAMtype BAM SortedByCoordinate \
			--outReadsUnmapped Fastx \
			--outFileNamePrefix ${sampleID}_${Lane}_
		

		fi
			



		"""
	
	}

	
