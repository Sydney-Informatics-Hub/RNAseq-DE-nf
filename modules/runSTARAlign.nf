process runSTARAlign {

        // where to publish the outputs
        tag "$sampleID runSTARAlign"
        publishDir "${params.outDir}/${sampleID}/STAR", mode:'copy'

input:

    path STAR_ref_index_path
    tuple val(sampleID), val(lane), val(runType), val(platform), val(sequencingCentre), val(library), path(r1Path), path(r2Path)

output:

    tuple val(sampleID), path ("${sampleID}_${lane}_Aligned.sortedByCoord.out.bam"), path ("${sampleID}_${lane}_SJ.out.tab") , emit: sample_lane_bam

script:

	def sampleID = sampleID
	def lane = lane
	def runType = runType
	def platform = platform
	def seqcentre = sequencingCentre
	def library = library
	def R1Path = r1Path ?: '' // Set default value for null R1 path
    def R2Path = r2Path ?: '' // Set default value for null R2 path
	

	"""
		# Mapping
		if [ ${runType} == 'PAIRED' ]; then

		 STAR \
            --runThreadN ${task.cpus} \
            --outBAMsortingThreadN ${task.cpus} \
            --genomeDir STARGeneratedIndexPath \
            --outBAMsortingBinsN 100 \
            --quantMode GeneCounts \
            --readFilesCommand zcat \
            --readFilesIn ${R1Path} ${R2Path} \
            --outSAMattrRGline ID:flowcell.${lane} PU:flowcell.${lane}.${sampleID} SM:${sampleID} PL:${platform} CN:${seqcentre} LB:${library} \
            --outSAMtype BAM SortedByCoordinate \
            --outReadsUnmapped Fastx \
            --outSAMunmapped Within KeepPairs \
            --outFileNamePrefix ${sampleID}_${lane}_

	else

	 STAR \
		--runThreadN ${task.cpus} \
		--outBAMsortingThreadN ${task.cpus} \
   		--genomeDir STARGeneratedIndexPath \
		--quantMode GeneCounts \
    	--outBAMsortingBinsN 100 \
		--readFilesCommand zcat \
		--readFilesIn ${R1Path} ${R2Path} \
		--outSAMattrRGline ID:flowcell.${lane} PU:flowcell.${lane}.${sampleID} SM:${sampleID} PL:${platform} CN:${seqcentre} LB:${library} \
		--outSAMtype BAM SortedByCoordinate \
		--outReadsUnmapped Fastx \
		--outFileNamePrefix ${sampleID}_${lane}_
		
	fi


	#touch ${sampleID}_${lane}_Aligned.sortedByCoord.out.bam
	#touch ${sampleID}_${lane}_SJ.out.tab

		"""
	
	}

	