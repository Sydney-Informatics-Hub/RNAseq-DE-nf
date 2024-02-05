process getMappingMetricRSeQC {

        // where to publish the outputs
        tag "$sampleID runSTARAlign"
        publishDir "${params.outDir}/STAR/${sampleID}/RSeQC", mode:'copy'

        input:
        	tuple val(sampleID) , file(R1_fastq), file(R2_fastq), val(seqcentre), val(platform), val(RUN_TYPE_SINGLE_PAIRED), val(lane),val(library)
		path refBed
		path sampleIDBam
		path sampleIDBamBai


	output:
		path ("${sampleID}_${lane}.infer_experiment.txt")
		//path ("${sampleID}_${lane}.read_distribution.txt")
		path ("${sampleID}_${lane}.bam_stat.txt")

        script:
	
		"""

		infer_experiment.py \
       			 -i $sampleIDBam \
        		 -r $refBed > ${sampleID}_${lane}.infer_experiment.txt

		#read_distribution.py \
		#	-i $sampleIDBam \
                #        -r $refBed > ${sampleID}_${lane}.read_distribution.txt

		bam_stat.py \
			-i $sampleIDBam > ${sampleID}_${lane}.bam_stat.txt

		


		"""
	
	}

	
