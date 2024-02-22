process convertGtfToBED {
  
        // where to publish the outputs
        tag "convertGtfToBED"
        publishDir "${params.outDir}", mode:'copy'

        input:
                path refGtf

        output:
                path ("Reference.bed")


        shell:
        """
        gtf2bed_SIH.pl ${refGtf} > "Reference.bed"

        """


}