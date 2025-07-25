nextflow.enable.dsl=2


process ExtractReads {
    errorStrategy 'ignore'
    
    publishDir params.output_gargammel, mode: 'symlink'

    input:
    val bam_file

    output:
    path ("*.fq.gz"), emit: read_file

    script:

    def output_prefix = "Observation_${bam_file.simpleName.replaceAll('_merged\\.fastq\\.bam$', '')}"

    """
    samtools1_21="/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/My_conda/bin/samtools"
    \${samtools1_21} fasta -@ {params.threads} -F 4  ${bam_file} | gzip >  "${output_prefix}.fa.gz" 
    """
}

process indexGenomesBowtie2 {
    // Define the input to be each genome
    
    input:
    val genome 

    // Specify the output files (index files will be generated by Bowtie2)
    
    output:
    path "*.rev.1.bt2" // This will output the Bowtie2 index files

    // Define the script to execute
    script:
    """
    # Run Bowtie2-build to index the genome
    ln -s ${params.ref_genomes}${genome}.fna . 
    bowtie2-build ${genome}.fna ${genome} &>> indexGenomes.log
    """
}

process indexGenomesBWA {
    // Define the input to be each genome
    input:
    val genome 

    // Specify the output files (index files will be generated by BWA)
    output:
    path "*.fna.sa" // This will output the BWA index files

    // Define the script to execute
    script:
    """
    # Run BWA to index the genome
    ln -s ${params.ref_genomes}${genome}.fna . 
    bwa index ${genome}.fna &>> indexGenomes.log
    """
}

process divergenceRunBowtie2 {

    publishDir params.output_gargammel, mode: 'symlink'  

    errorStrategy 'ignore'

    conda "${params.conda_env2}"

    maxForks 499

    input:
    tuple val(fq_file), val(Indexed_genomes), val(sizes)

    output:
    file "*.bam"
    file "*.pdf"
    file "*.txt"
    file "*.amber"
    file "*.png"
    file "*.csv"

    script:
    def base_name_only = fq_file.simpleName.replaceAll('\\.fa\\.gz$', '') // Remove the suffix to get the base name

    // Run Bowtie2
    """
    echo "Running divergence script and Bowtie2 alignment for ${fq_file}..."

    samtools1_21="/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/My_conda/bin/samtools"

    base_name=\$(basename ${fq_file} .fa.gz)
    ref_species=\$(basename ${Indexed_genomes})

    # Map the reads to the reference using Bowtie2
    bowtie2 --quiet --very-sensitive --threads ${params.threads} -x ${Indexed_genomes}  -U ${fq_file} | \
    \${samtools1_21} view -b -F 4 -@ ${params.threads} | \
    \${samtools1_21} sort -@ ${params.threads} -O bam -o "\${base_name}_ref-\${ref_species}_ProgBowtie2.bam" &>> divergenceRun.log

    # Subsample the desired number of reads for Bowtie2 output
    output_sampled_bam="\${base_name}_ref-\${ref_species}_nread${sizes}.bam"
    \${samtools1_21} view -H "\${base_name}_ref-\${ref_species}_ProgBowtie2.bam" > header.sam   # Extract the header
    \${samtools1_21} view "\${base_name}_ref-\${ref_species}_ProgBowtie2.bam" | shuf -n ${sizes} > random_reads.sam  # Extract N random reads
    cat header.sam random_reads.sam | \${samtools1_21} view -bS - > \${output_sampled_bam}   # Combine header and random reads into a new BAM file
    \${samtools1_21} sort \${output_sampled_bam} -o "\${output_sampled_bam%.bam}_ProgBowtie2_sorted.bam" &>> divergenceRun.log

    # Index the resulting BAM file
    \${samtools1_21} index -@ ${params.threads} "\${output_sampled_bam%.bam}_ProgBowtie2_sorted.bam" &>> divergenceRun.log

    # Run Amber for Bowtie2 output
    python3 /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/Amber_last_noplot.py -b "\${output_sampled_bam%.bam}_ProgBowtie2_sorted.bam" -o "\${output_sampled_bam%.bam}_ProgBowtie2_sorted.amber" &>> divergenceRun.log
    java -jar /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/DamageProfiler-1.1-java11.jar -i "\${output_sampled_bam%.bam}_ProgBowtie2_sorted.bam" -o DamageProfiler_output &>> divergenceRun.log
    mv DamageProfiler_output/5p_freq_misincorporations.txt "./\${base_name}_ref-\${ref_species}_nread${sizes}_ProgBowtie2_5p_freq_misincorporations.txt" &>> divergenceRun.log
    mv DamageProfiler_output/3p_freq_misincorporations.txt "./\${base_name}_ref-\${ref_species}_nread${sizes}_ProgBowtie2_3p_freq_misincorporations.txt" &>> divergenceRun.log
    mv DamageProfiler_output/DamagePlot.pdf "./\${base_name}_ref-\${ref_species}_nread${sizes}_ProgBowtie2_DamagePlot.pdf" &>> divergenceRun.log
    rm -r DamageProfiler_output &>> divergenceRun.log
    rm *sorted.bam.bai *random_reads.sam *header.sam &>> divergenceRun.log
    rm "\${base_name}_ref-\${ref_species}_ProgBowtie2.bam" &>> divergenceRun.log

    # Run Pydamage for Bowtie2 output
    \${samtools1_21} index "\${output_sampled_bam%.bam}_ProgBowtie2_sorted.bam"
    pydamage analyze "\${output_sampled_bam%.bam}_ProgBowtie2_sorted.bam" -p ${params.threads} -pl -g 
    cp pydamage_results/pydamage_results.csv . 
    cp pydamage_results/plots/reference.png . 
    mv pydamage_results.csv "Pydamage_\${base_name}_ref-\${ref_species}_nread${sizes}_ProgBowtie2.csv"
    mv reference.png "Pydamage_\${base_name}_ref-\${ref_species}_nread${sizes}_ProgBowtie2.png"
    rm -r pydamage_results
    python3 /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/Plot_identity_and_damage.py "\${output_sampled_bam%.bam}_ProgBowtie2_sorted.bam" "./\${base_name}_ref-\${ref_species}_nread${sizes}_ProgBowtie2_3p_freq_misincorporations.txt" "./\${base_name}_ref-\${ref_species}_nread${sizes}_ProgBowtie2_5p_freq_misincorporations.txt" "\${output_sampled_bam%.bam}_ProgBowtie2_sorted.amber" "Pydamage_\${base_name}_ref-\${ref_species}_nread${sizes}_ProgBowtie2.csv" &>> divergenceRun.log 
    """
}

process divergenceRunBWA {

    errorStrategy 'ignore'

    publishDir params.output_gargammel, mode: 'symlink'  

    conda "${params.conda_env2}"

    maxForks 499

    input:
    tuple val(fq_file), val(Indexed_genomes), val(sizes)

    output:
    file "*.bam"
    file "*.pdf"
    file "*.txt"
    file "*.amber"
    file "*.png"
    file "*.csv"

    script:
    def base_name_only = fq_file.simpleName.replaceAll('\\.fa\\.gz$', '') // Remove the suffix to get the base name

    """
    echo "Running divergence script and BWA alignment for ${fq_file}..."

    samtools1_21="/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/My_conda/bin/samtools"

    base_name=\$(basename ${fq_file} .fa.gz)
    ref_species=\$(basename ${Indexed_genomes})

    # Map the reads to the reference using BWA
    bwa aln -l 16500 -n 0.01 -o 2 -t ${params.threads} ${Indexed_genomes}.fna ${fq_file} > "\${base_name}_ref-\${ref_species}_ProgBWA.sai"
    bwa samse ${Indexed_genomes}.fna "\${base_name}_ref-\${ref_species}_ProgBWA.sai" ${fq_file} | \
    \${samtools1_21} view -b -F 4 -@ ${params.threads} | \
    \${samtools1_21} sort -@ ${params.threads} -O bam -o "\${base_name}_ref-\${ref_species}_ProgBWA.bam" &>> divergenceRun.log

    # Subsample the desired number of reads for BWA output
    output_sampled_bam_bwa="\${base_name}_ref-\${ref_species}_nread${sizes}.bam"
    \${samtools1_21} view -H "\${base_name}_ref-\${ref_species}_ProgBWA.bam" > header.sam   # Extract the header
    \${samtools1_21} view "\${base_name}_ref-\${ref_species}_ProgBWA.bam" | shuf -n ${sizes} > random_reads.sam  # Extract N random reads
    cat header.sam random_reads.sam | \${samtools1_21} view -bS - > \${output_sampled_bam_bwa}   # Combine header and random reads into a new BAM file
    \${samtools1_21} sort \${output_sampled_bam_bwa} -o "\${output_sampled_bam_bwa%.bam}_ProgBWA_sorted.bam" &>> divergenceRun.log

    # Index the resulting BAM file
    \${samtools1_21} index -@ ${params.threads} "\${output_sampled_bam_bwa%.bam}_ProgBWA_sorted.bam" &>> divergenceRun.log

    # Run Amber for BWA output
    python3 /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/Amber_last_noplot.py -b "\${output_sampled_bam_bwa%.bam}_ProgBWA_sorted.bam" -o "\${output_sampled_bam_bwa%.bam}_ProgBWA_sorted.amber" &>> divergenceRun.log
    java -jar /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/DamageProfiler-1.1-java11.jar -i "\${output_sampled_bam_bwa%.bam}_ProgBWA_sorted.bam" -o DamageProfiler_output &>> divergenceRun.log
    mv DamageProfiler_output/5p_freq_misincorporations.txt "./\${base_name}_ref-\${ref_species}_nread${sizes}_ProgBWA_5p_freq_misincorporations.txt" &>> divergenceRun.log
    mv DamageProfiler_output/3p_freq_misincorporations.txt "./\${base_name}_ref-\${ref_species}_nread${sizes}_ProgBWA_3p_freq_misincorporations.txt" &>> divergenceRun.log
    mv DamageProfiler_output/DamagePlot.pdf "./\${base_name}_ref-\${ref_species}_nread${sizes}_ProgBWA_DamagePlot.pdf" &>> divergenceRun.log
    rm -r DamageProfiler_output &>> divergenceRun.log
    rm *sorted.bam.bai *random_reads.sam *header.sam &>> divergenceRun.log
    rm "\${base_name}_ref-\${ref_species}_ProgBWA.bam" &>> divergenceRun.log

    # Run Pydamage for BWA output
    \${samtools1_21} index "\${output_sampled_bam_bwa%.bam}_ProgBWA_sorted.bam"
    pydamage analyze "\${output_sampled_bam_bwa%.bam}_ProgBWA_sorted.bam" -p ${params.threads} -pl -g 
    cp pydamage_results/pydamage_results.csv . 
    cp pydamage_results/plots/reference.png . 
    mv pydamage_results.csv "Pydamage_\${base_name}_ref-\${ref_species}_nread${sizes}_ProgBWA.csv"
    mv reference.png "Pydamage_\${base_name}_ref-\${ref_species}_nread${sizes}_ProgBWA.png"
    rm -r pydamage_results
    python3 /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/Plot_identity_and_damage.py "\${output_sampled_bam_bwa%.bam}_ProgBWA_sorted.bam" "./\${base_name}_ref-\${ref_species}_nread${sizes}_ProgBWA_3p_freq_misincorporations.txt" "./\${base_name}_ref-\${ref_species}_nread${sizes}_ProgBWA_5p_freq_misincorporations.txt" "\${output_sampled_bam_bwa%.bam}_ProgBWA_sorted.amber" "Pydamage_\${base_name}_ref-\${ref_species}_nread${sizes}_ProgBWA.csv" &>> divergenceRun.log
    """
}


workflow {

    // Create a channel from BAM files containing "udg"

    fq_file = Channel.fromPath("${params.output_fastq}/${params.species}*merged.fastq.gz")
 
    // Create a channel from the genome list
    Genome_list = Channel.from(params.genome_list)
    
    // Index the genomes
    indexGenomesBowtie2(Genome_list)
    indexGenomesBWA(Genome_list)
    
    //Indexed_genomes = indexGenomes.out.flatten().collect { it.toString().split('\\.')[0] }
    //indexGenomesBowtie2.out.view()
    //indexGenomesBWA.out.view()
    Indexed_genomesBowtie2 = indexGenomesBowtie2.out.flatten().map { it.toString().split('\\.')[0] }
    Indexed_genomesBWA = indexGenomesBWA.out.flatten().map { "${it.toString().split('\\.')[0]}" }
    //Indexed_genomesBowtie2.view { "Bowtie2: $it" }
    //Indexed_genomesBWA.view { "BWA: $it" }

    sizes = Channel.of(10000,5000, 1000, 500, 100)
    //sizes = Channel.of(10000,100)

    divergence_data_Bowtie2 = fq_file.combine(Indexed_genomesBowtie2).combine(sizes)

    divergence_data_BWA = fq_file.combine(Indexed_genomesBWA).combine(sizes)

    divergence_data_BWA = divergence_data_BWA.mix(divergence_data_BWA)
    divergence_data_Bowtie2 = divergence_data_Bowtie2.mix(divergence_data_Bowtie2)

    // Run divergenceRun for output from ExtractReads and gargammelRun
    divergenceRunBowtie2(divergence_data_Bowtie2)
    //divergenceRunBWA(divergence_data_BWA)
}
