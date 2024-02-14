process makeSTARIndex {
  
        // where to publish the outputs
        tag "makeSTARIndex"
        publishDir "${params.outDir}/INDEX/STAR", mode:'copy'

        input:
        	path referenceFasta
		path refGtf
		val(NCPUS)		

        output:
           path ("STARGeneratedIndexPath")

        script:

        """
		
	STAR \
		--runThreadN ${NCPUS} \
        	--runMode genomeGenerate \
       		--genomeDir STARGeneratedIndexPath \
        	--genomeFastaFiles ${referenceFasta} \
        	--sjdbGTFfile ${refGtf} \
        	--sjdbOverhang 99

	"""




	}


