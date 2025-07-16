// Create a channel from the input file (assuming it's a TSV file)
Channel
    .fromPath(params.inputFile)               // Load the input file path from params
    .splitCsv(sep: '\t', header: true)        // Split it into rows as a TSV file, setting header as true
    .map { row -> [accession: row.archive_data_accession, treatment: row.library_treatment] }  // Extract the 'archive_data_accession' column
    .set { accessionsAndTreatments }                      // Save the channel as accessionsAndTreatments


accessionsAndTreatments.view()

process DOWNLOAD_FASTQ {

    //errorStrategy 'ignore'

    conda '/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/My_conda'

    publishDir params.output_reads, mode: 'symlink'

    input:
    tuple val(accession), val(treatment) 

    output:
    path("*_merged.fastq.gz"), emit: merged_reads 

    script:
    """
    echo ""
    echo "############################################################################################################"
    echo "###    Downloading reads of : ${accession} ${treatment}                                               "
    echo "###    All output will be written in :  ${params.output_reads}                                             "
    /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/sratoolkit.3.1.0-ubuntu64/bin/fasterq-dump --threads ${params.threads_DOWNLOAD_FASTQ} --outdir ./ ${accession} &>> DOWNLOAD_FASTQ.log
    pigz -p ${params.threads_DOWNLOAD_FASTQ} ./${accession}_1.fastq &>> DOWNLOAD_FASTQ.log
    pigz -p ${params.threads_DOWNLOAD_FASTQ} ./${accession}_2.fastq &>> DOWNLOAD_FASTQ.log
    for file in *.gz; do mv "\${file}" "${params.sp_name}_\${file}"; done
    echo "Files downloaded successfully"
    echo ""
    echo "### Filtering and merging reads"
    echo "Filtering reads of : ${accession}"
    filename=${params.sp_name}_\$(basename ${accession} _1.fastq.gz)
    time_start=\$(date +%s)
    mkdir -p ${params.output_reads}
    echo ""
    fastp -i \${filename}_1.fastq.gz  -I \${filename}_2.fastq.gz -m --merged_out \${filename}_merged.fastq.gz --dedup -m 3 -o \${filename}_trimmed.R1.fastq.gz -O \${filename}_trimmed.R2.fastq.gz -h \${filename}_fastp.html -l 30 -w ${params.threads_DOWNLOAD_FASTQ} &>> DOWNLOAD_FASTQ.log
    echo "Reads filtered successfully"
    mv "\${filename}_merged.fastq.gz" "${params.sp_name}_${accession}_${treatment}_merged.fastq.gz"
    #
    rm *_1.fastq.gz *_2.fastq.gz *_trimmed.R1.fastq.gz *_trimmed.R2.fastq.gz
    """
}


process MAP_READS {
    // Set the maximum number of threads to use for this process
    errorStrategy { task.exitStatus == 143 ? 'ignore' : 'terminate' }

    publishDir params.output_mappings, mode: 'symlink'    

    // 3.1 Define the input file
    input:
    each (merged_reads)

    // 3.2 Define the output file
    output:
    path("*.bam"), emit: mapped_reads

    // 3.3 Define the script
    script:
    """
    samtools1_21="/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/My_conda/bin/samtools"
    ref_assembly="${params.dir_ref}/${params.sp_name}.fna"
    merged_reads_name=\$(basename ${merged_reads.baseName} _merged.fastq.gz)
    #merged_reads="${params.dir_bam}/ALL_merged.Yersina_control.fastq.gz"
    echo indexing
    bowtie2-build --quiet --large-index --threads ${params.threads_MAP_READS} -f \${ref_assembly} ${params.sp_name} &>> MAP_READS.log
    echo mapping
    bowtie2 --quiet --very-sensitive --threads ${params.threads_MAP_READS} -x ${params.sp_name} -U ${merged_reads} | \${samtools1_21} view --verbosity 0 -b -F 4 -@ ${params.threads_MAP_READS} | \${samtools1_21} sort --verbosity 0 -@ ${params.threads_MAP_READS} -O bam -o \${merged_reads_name}.bam  &>> MAP_READS.log
    echo indexing 
    \${samtools1_21} index -@ ${params.threads_MAP_READS} \${merged_reads_name}.bam  &>> MAP_READS.logy
    """
}



// 4. Run the workflow and connect the two processes
workflow {

    // Read accessions from the input file (already done before the workflow block)
    
    DOWNLOAD_FASTQ(accessionsAndTreatments)

    merged_reads=DOWNLOAD_FASTQ.out.merged_reads.flatten()

    MAP_READS(merged_reads)

}
