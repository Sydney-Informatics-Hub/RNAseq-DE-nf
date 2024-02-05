




	input:



	output:


	


	shell:

	'''

	# Read the sample sheet CSV file and extract unique sample IDs
	sampleIDs=$(awk -F ',' 'NR>1 {print $1}' /scratch/er01/ndes8648/pipeline_work/nextflow/INFRA-121-RNASeq-DE/RNAseq-DE-nf/sampleSheet_both_single_paired.csv | sort -u)


	# Initialize an empty array to store the sample IDs
	declare -a sampleID_array

	# Loop through each sample ID and add it to the array
	while IFS= read -r line; do
    		sampleID_array+=("$line")
	done <<< "$sampleIDs"

	# Print the array
	echo "Sample IDs:"
	echo "(""${sampleID_array[*]}"")"



	'''


	
