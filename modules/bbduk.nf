// Define the process
process bbduk {	

    debug = true //turn to false to stop printing command stdout to screen
    publishDir "${params.outDir}/${sampleID}/bbduk/trimmed", mode: 'copy'

    input: 
	path adapters_bbmap
    	tuple val(sampleID), val(Lane), path(R1), path(R2), val(SEQUENCING_CENTR), val(PLATFORM), val(RUN_TYPE_SINGLE_PAIRED), val(LIBRARY)
    	val(NCPUS)

    output:
    	path("${sampleID}.R1.trimmed.fastq.gz")
    	path("${sampleID}.R2.trimmed.fastq.gz") , optional: true

    script:
    """
	readlen=100
	
	if [ "${RUN_TYPE_SINGLE_PAIRED}" == 'PAIRED' ]; then
	
	bbduk.sh -Xmx6g \
		threads=${NCPUS} \
		in=${R1} \
		in2=${R2} \
		ref=${adapters_bbmap} \
		out=${sampleID}.R1.trimmed.fastq.gz \
    		out2=${sampleID}.R2.trimmed.fastq.gz \
		ktrim=r \
		k=23 \
		mink=11 \
		hdist=1 \
		tpe \
		tbo \
		overwrite=true \
		trimpolya=readlen

	else
	
	bbduk.sh -Xmx6g \
                threads=${NCPUS} \
                in=${R1} \
                ref=${adapters_bbmap} \
                out=${sampleID}.R1.trimmed.fastq.gz \
                ktrim=r \
                k=23 \
                mink=11 \
                hdist=1 \
                tpe \
                tbo \
                overwrite=true \
                trimpolya=readlen	

	fi	
	
    """
}
