// Create a channel from the input file (assuming it's a TSV file)
Channel
    .fromPath(params.inputFile)               // Load the input file path from params
    .splitCsv(sep: '\t', header: true)        // Split it into rows as a TSV file, setting header as true
    .map { row -> [accession: row.archive_data_accession, treatment: row.library_treatment, library: row.library_layout, link: row.download_links, archive: row.archive] }  // Extract the 'archive_data_accession' column
    .set { accessionsAndTreatments }                      // Save the channel as accessionsAndTreatments


accessionsAndTreatments.view()

process DOWNLOAD_FASTQ {

    publishDir params.output_reads, mode: 'symlink'

    input:
    tuple val(accession), val(treatment), val(layout), val(link), val(archive)

    output:
    path("*merged.fastq.gz"), emit: merged_reads 

    script:
    """
    echo ""
    echo "############################################################################################################"
    echo "###    Downloading reads of : ${accession} ${treatment}                                               "
    echo "###    All output will be written in :  ${params.output_reads}                                             "

    # Download reads based on layout
    #/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/sratoolkit.3.1.0-ubuntu64/bin/fasterq-dump --threads ${params.threads_DOWNLOAD_FASTQ} --outdir ./ ${accession} &>> Download_fast.log

    if [[ "${layout}" == "PAIRED" ]]; then
        /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/sratoolkit.3.1.0-ubuntu64/bin/fasterq-dump --threads ${params.threads_DOWNLOAD_FASTQ} --outdir ./ ${accession} &>> Download_fast.log
        pigz -p ${params.threads_DOWNLOAD_FASTQ} ./${accession}_1.fastq &>> Download_fast.log
        pigz -p ${params.threads_DOWNLOAD_FASTQ} ./${accession}_2.fastq &>> Download_fast.log
        for file in *.gz; do mv "\${file}" "${params.sp_name}_\${file}"; done
        filename=${params.sp_name}_\$(basename ${accession} _1.fastq.gz)

        echo "Files downloaded successfully"
        echo ""
        echo "### Filtering and merging paired-end reads"
        fastp -i \${filename}_1.fastq.gz  -I \${filename}_2.fastq.gz -m --merged_out \${filename}_merged.fastq.gz --dedup -m 3 -o \${filename}_trimmed.R1.fastq.gz -O \${filename}_trimmed.R2.fastq.gz -h \${filename}_fastp.html -l 30 -w ${params.threads_DOWNLOAD_FASTQ}  &>> Download_fast.log
        rm *_1.fastq.gz *_2.fastq.gz *_trimmed.R1.fastq.gz *_trimmed.R2.fastq.gz 
    else
        if [[ "${archive}" == "ENA" ]]; then
           /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/sratoolkit.3.1.0-ubuntu64/bin/fasterq-dump --threads ${params.threads_DOWNLOAD_FASTQ} --outdir ./ ${accession} &>> Download_fast.log
           pigz -p ${params.threads_DOWNLOAD_FASTQ} ./${accession}.fastq &>> Download_fast.log
        else 
           wget ${link}
        fi 

        mv "${accession}.fastq.gz" "${params.sp_name}_${accession}_${treatment}_single.fastq.gz"
        filename=${params.sp_name}_${accession}_${treatment}_single

        echo "Files downloaded successfully"
        echo ""
        echo "### Filtering single-end reads"
        fastp -i \${filename}.fastq.gz -o \${filename}_merged.fastq.gz  -l 30 -w ${params.threads_DOWNLOAD_FASTQ} &>> Download_fast.log
        rm \${filename}.fastq.gz
    fi

    echo "Reads filtered successfully"
    mv "\${filename}_merged.fastq.gz" "${params.sp_name}_${accession}_${treatment}_merged.fastq.gz"

    """
}

 
// 4. Run the workflow and connect the two processes
workflow {

    // Read accessions from the input file (already done before the workflow block)
    
    DOWNLOAD_FASTQ(accessionsAndTreatments)

}
