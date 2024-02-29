process createSalmonPCA {
  
        // where to publish the outputs
        tag "createSalmonPCA"
        publishDir "${params.outDir}", mode:'copy'

        input:
                path merged_salmon_counts
                path samples_info


        output:
                path ("pca_plot_salmon.png")

        script:
        """
        salmon_PCA.r ${merged_salmon_counts} ${samples_info}

        """

}