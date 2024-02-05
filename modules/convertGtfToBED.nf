process convertGtfToBED {
  
        // where to publish the outputs
        tag "convertGtfToBED"
        publishDir "${params.outDir}", mode:'copy'

        input:
                path refGtf


        output:
                path ("Reference.bed")








	shell:
	'''

	# Input GTF file
	input_gtf=!{refGtf}

	# Output BED file
	output_bed="Reference.bed"

	# Convert GTF to BED
	awk -F '\t' 'BEGIN {OFS="\t"} $3=="exon" {print $1,$4-1,$5,$9,0,$7}' "$input_gtf" > "$output_bed"

	#echo "Conversion complete. BED file saved as: $output_bed"

	'''



	}
