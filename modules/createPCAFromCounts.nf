process createPCAFromCounts {
  
        // where to publish the outputs
        tag "createPCAFromCounts"
        publishDir "${params.outDir}", mode:'copy'

        input:
                path merged_counts_STAR
                path samples_info


        output:
                path ("pca_plot_STAR.png")

        script:
        """
        create_PCA.r ${merged_counts_STAR} ${samples_info}

        """

}