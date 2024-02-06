process getMappingMetricRSeQC {

        // where to publish the outputs
        tag "$uniqueSampleID getMappingMetricRSeQC"
        publishDir "${params.outDir}/STAR/${uniqueSampleID}/RSeQC", mode:'copy'

        input:
		val uniqueSampleID
		path refBed
		path uniqueSampleIDBam
		path uniqueSampleIDBamBai



	output:
		val ("$uniqueSampleID")
		path ("${uniqueSampleID}.infer_experiment.txt")
		//path ("${uniqueSampleID}.read_distribution.txt")
		path ("${uniqueSampleID}.bam_stat.txt")
			


        script:
	
		"""

		infer_experiment.py \
       			 -i $uniqueSampleIDBam \
        		 -r $refBed > ${uniqueSampleID}.infer_experiment.txt

		#read_distribution.py \
		#	-i $uniqueSampleIDBam \
                #        -r $refBed > ${uniqueSampleID}.read_distribution.txt

		bam_stat.py \
			-i $uniqueSampleIDBam > ${uniqueSampleID}.bam_stat.txt

		


		"""
	
	}

	
