process makeSTARIndex {
  
        // where to publish the outputs
        tag "makeSTARIndex ${params.refFasta}"
        publishDir "${params.outDir}", mode:'copy'

        input:
        path referenceFasta
		path refGtf	

        output:
        path ("STARGeneratedIndexPath"), emit: STAR_ref_index_path

        script:

		def refFasta = referenceFasta
        """
		STAR \
			--runThreadN ${task.cpus} \
        	--runMode genomeGenerate \
       		--genomeDir STARGeneratedIndexPath \
        	--genomeFastaFiles ${refFasta} \
        	--sjdbGTFfile ${refGtf} \
        	--sjdbOverhang 99
	"""
	}