// Define the process
process bbduk {	

    debug = true //turn to false to stop printing command stdout to screen
    publishDir "${params.outDir}/${sampleID}/bbduk/trimmed", mode: 'copy'

    input: 
	path adapters_bbmap
    	tuple val(sampleID), val(Lane), path(R1), path(R2), val(SEQUENCING_CENTR), val(PLATFORM), val(RUN_TYPE_SINGLE_PAIRED), val(LIBRARY)

    output:
    	tuple val(sampleID), val(Lane), val(RUN_TYPE_SINGLE_PAIRED), val(PLATFORM), val(SEQUENCING_CENTR), val(LIBRARY), path("${sampleID}_${Lane}.R1.trimmed.fastq.gz"), path("${sampleID}_${Lane}.R2.trimmed.fastq.gz"), emit: trimmed_fq

    script:
    """
	readlen=100
	
	if [ "${RUN_TYPE_SINGLE_PAIRED}" == 'PAIRED' ]; then
	
	bbduk.sh -Xmx6g \
		threads=${task.cpus} \
		in=${R1} \
		in2=${R2} \
		ref=${adapters_bbmap} \
		out=${sampleID}_${Lane}.R1.trimmed.fastq.gz \
    	out2=${sampleID}_${Lane}.R2.trimmed.fastq.gz \
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
                threads=${task.cpus} \
                in=${R1} \
                ref=${adapters_bbmap} \
                out=${sampleID}_${Lane}.R1.trimmed.fastq.gz \
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