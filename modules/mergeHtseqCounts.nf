process mergeHtseqCounts {
  
        // where to publish the outputs
        tag "$sampleID runSTARAlign"
        publishDir "${params.outDir}", mode:'copy'

        input:
                path(allHtseqCountFiles)


        output:
                path ("merged_counts_STAR.txt")

        shell:

        '''
	# Output file name
        output_file="merged_counts_STAR.txt"

        # Initialize an associative array to store counts
        declare -A counts


	# Collect sample IDs
        sample_ids=()
        for sample_dir in !{allHtseqCountFiles}; do
        #if [[ -d "$sample_dir" ]]; then
                sample_id=$(basename "$sample_dir")
                sample_ids+=("$sample_id")
        #fi
        done

	# Loop through each .counts file
	for counts_file in !{allHtseqCountFiles}; do
    		# Get sample ID from counts file path
    		#sample_id=$(basename "$(dirname "$counts_file")")
		sample_id=$(basename "$sample_dir")

    		# Read each line of the file
    		while IFS=$'\t' read -r gene_id count; do
        		# Add count to existing gene_id or initialize it
        		counts["$gene_id"]+="\t$count"
    		done < "$counts_file"
	done



	# Write merged counts to output file
	#echo -e "GENEID\t${sample_ids[*]}" > "$output_file"
	echo -e "GENEID\t${sample_ids[*]// /\\t}" > "$output_file"
	for gene_id in "${!counts[@]}"; do
		    echo -e "$gene_id${counts[$gene_id]}" >> "$output_file"
	done


	'''

        }


