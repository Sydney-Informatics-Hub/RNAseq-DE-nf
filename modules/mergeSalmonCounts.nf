process mergeSalmonCounts {
  
        debug = true //turn to false to stop printing command stdout to screen
    publishDir "${params.outDir}", mode: 'copy'   

        input:
            path(counts_gene)


        output:
                path("final_merged_counts_Salmon.tsv"), emit: merged_salmon_counts

        shell:

    '''
   
#!/usr/bin/bash

# Output file name
output_file="merged_counts_Salmon.tsv"

# Collect sample IDs
sample_ids=()
for sample_dir in !{counts_gene}; do
    #if [[ -d "$sample_dir" ]]; then
    sample_id=$(basename "$sample_dir")
    sample_ids+=("$sample_id")
    #fi
done

# Print custom header
echo -e "GENEID\t${sample_ids[*]}" > "$output_file"

# Merge counts and append to output file
awk 'NR == 1 { next } FNR > 1 { a[$1] = a[$1]"\t"$3 } END { for (i in a) print i a[i] }' !{counts_gene} >> "$output_file"

# Confirm the space between columns
tr -s ' ' '\t' < ${output_file} > final_${output_file}
    
    '''
}