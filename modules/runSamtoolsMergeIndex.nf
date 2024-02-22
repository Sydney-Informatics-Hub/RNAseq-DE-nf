process runSamtoolsMergeIndex {


	// where to publish the outputs
        tag "runSamtoolsMergeIndex"
        publishDir "${params.outDir}/${sampleID}/STAR", mode:'copy'

	input:
        tuple val(sampleID), path(path_bam), path(path_SJ)
		

	output:
	    tuple val(sampleID), path ("${sampleID}_final.bam") , path ("${sampleID}_final.bam.bai") , emit: final_bam
		

	script:
	def sampleID = sampleID
 	def path_bam = path_bam
  	def path_SJ = path_SJ

	"""
	 samtools merge -f -@ ${task.cpus} ${sampleID}_final.bam ${path_bam}

	 samtools index -@ ${task.cpus} ${sampleID}_final.bam -o ${sampleID}_final.bam.bai

	#touch ${sampleID}_final.bam
	#touch ${sampleID}_final.bam.bai

	"""

	

}