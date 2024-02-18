process getMappingMetricRSeQC {
  
        // where to publish the outputs
        tag "$sampleID getMappingMetricRSeQC"
        publishDir "${params.outDir}/${sampleID}/STAR/RSeQC", mode:'copy'

        input:
                tuple val(sampleID), path(final_bam_file), path(final_bam_index)
                path refBed

        output:
                tuple val(sampleID), path ("${sampleID}.infer_experiment.txt"), path ("${sampleID}.read_distribution.txt"), path ("${sampleID}.bam_stat.txt"), emit: star_metrics
        
        script:

                """

                infer_experiment.py \
                         -i $final_bam_file \
                         -r $refBed > ${sampleID}.infer_experiment.txt

                read_distribution.py \
                       -i $final_bam_file \
                        -r $refBed > ${sampleID}.read_distribution.txt

                bam_stat.py \
                        -i $final_bam_file > ${sampleID}.bam_stat.txt

                """

        }

