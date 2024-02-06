process runSamtoolsMergeIndex {


	// where to publish the outputs
        tag "$uniqueSampleID runSamtoolsMergeIndex"
        publishDir "${params.outDir}/STAR/${uniqueSampleID}", mode:'copy'

	input:
                each uniqueSampleID
                path(sampleID_bams)
                val(NCPUS)

	output:
		path ("${uniqueSampleID}.final.bam")
		path ("${uniqueSampleID}.final.bam.bai")
		val ("$uniqueSampleID")

	shell:
	'''

	#https://github.com/nextflow-io/nextflow/issues/3962

	sampbams=()
	sampbams+=("$(ls !{uniqueSampleID}*_Aligned.sortedByCoord.out.bam)")	


	echo My working DIR: $PWD
	echo bams: $sampbams
	
	final=!{uniqueSampleID}.final.bam

	# Define an empty array
	sampbams_final=()

	# Iterate over each element in sampbams array
		for filename in \${sampbams[@]}; do
   		# Concatenate $PWD/ to each filename and add it to sampbams_final array
    			sampbams_final+=(\$PWD/$filename)
		done


	echo bams final: $sampbams_final


	samtools merge -f -@ !{NCPUS} "$final" "$sampbams_final"

	samtools index -@ !{NCPUS} "$final" -o "$final".bai


	'''

	

}
