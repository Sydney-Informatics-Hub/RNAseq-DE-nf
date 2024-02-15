process runSTARAlign {

        // where to publish the outputs
        tag "$sampleID runSTARAlign"
        publishDir "${params.outDir}/${sampleID}/STAR", mode:'copy'

input:
// INPUT ORDER HERE NEEDS TO MATCH DEF ORDER BELOW 
    path STAR_ref_index_path
    tuple val(sampleID), val(lane), val(runType), val(platform), val(sequencingCentre), val(library), path(r1Path), path(r2Path)

output:
/// I CHANGED THIS TO CAPTURE SINGLE AND PAIRED END READS 
/// AS IT WAS, IT FAILS WHEN R2 IS NULL AS YOU HAD SPECIFIED $LANE IN THE OUTPUT PATH
    path ("${sampleID}_*Aligned.sortedByCoord.out.bam") , emit: sampleID_lane_bam
	path ("${sampleID}_*SJ.out.tab")			   , emit: sampleID_lane_SJ_tab

script:
	/// ORDER OF THESE DEFINITIONS IS BASED ON INPUT ORDER ABOVE
	/// YOU CAN TEST THIS BY SWAPPING val(runType) WITH val(platform) IN THE INPUT SECTION
	/// AND SEEING HOW THE SCRIPT FAILS
	def sampleID = sampleID
	def Lane = lane
	def runType = runType
	def platform = platform
	def seqcentre = sequencingCentre
	def library = library
	def R1Path = r1Path ?: '' // Set default value for null R1 path
    def R2Path = r2Path ?: '' // Set default value for null R2 path
	
	/// NOTE IN MAIN.NF PROCESS DEFINITIONS AND SOME PROCESS FILES I WAS ADJUSTING 
	/// THAT I HAVE REMOVE ${NCPUS}. THERE IS NO NEED TO FILL THIS IN MANUALLY. 
	/// SEE: https://www.nextflow.io/docs/latest/process.html FOR HOW TO USE TASK.CPU, TASK.MEMORY ETC 
	/// TO CONNECT CONFIGURATION TO YOUR PROCESSES.
	"""
		# Mapping
		if [ ${runType} == 'PAIRED' ]; then

		STAR \
            --runThreadN ${task.cpus} \
            --outBAMsortingThreadN ${task.cpus} \
            --genomeDir STAR_INDEX/STARGeneratedIndexPath \
            --outBAMsortingBinsN 100 \
            --quantMode GeneCounts \
            --readFilesCommand zcat \
            --readFilesIn ${R1Path} ${R2Path} \
            --outSAMattrRGline ID:flowcell.${Lane} PU:flowcell.${Lane}.${sampleID} SM:${sampleID} PL:${platform} CN:${seqcentre} LB:${library} \
            --outSAMtype BAM SortedByCoordinate \
            --outReadsUnmapped Fastx \
            --outSAMunmapped Within KeepPairs \
            --outFileNamePrefix ${sampleID}_${Lane}_

	else

	STAR \
		--runThreadN ${task.cpus} \
		--outBAMsortingThreadN ${task.cpus} \
   		--genomeDir STAR_INDEX/STARGeneratedIndexPath \
		--quantMode GeneCounts \
    	--outBAMsortingBinsN 100 \
		--readFilesCommand zcat \
		--readFilesIn ${R1Path} ${R2Path} \
		--outSAMattrRGline ID:flowcell.${Lane} PU:flowcell.${Lane}.${sampleID} SM:${sampleID} PL:${platform} CN:${seqcentre} LB:${library} \
		--outSAMtype BAM SortedByCoordinate \
		--outReadsUnmapped Fastx \
		--outFileNamePrefix ${sampleID}_${Lane}_
		
	fi
		"""
	
	}

	