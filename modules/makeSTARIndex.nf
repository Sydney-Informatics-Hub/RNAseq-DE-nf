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
		/// NEED TO BE CONSISTENT IN NAMING BETWEEN INPUT CHANNEL AND VARIABLE DEFINITION
		/// IF YOU USE refFasta IN INPUT CHANNEL, YOU NEED TO USE refFasta IN VARIABLE DEFINITION
		/// OTHERWISE, IT IS NOT VISIBLE TO THE SCRIPT BLOCK
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