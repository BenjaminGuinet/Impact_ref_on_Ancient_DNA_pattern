nextflow.enable.dsl=2

process gargammelRun {

    errorStrategy 'ignore'

    publishDir params.output_gargammel, mode: 'symlink'  
   
    conda "${params.conda_env1}"

    maxForks 499

    input:
    tuple  val(damage_param), val(bam_file), val(read_length_param)

    output:
    path ("*.fa.gz"), emit: gargammel_read_file

    // Set the damage level in the Groovy context, before the script block
    when : 
        damage_level = 
        (damage_param == " ") ? "0DOperc" :
        //(damage_param == "-mf misincorporation.txt") ? "normDOperc" :
        (damage_param == "-m b,0.024,0.36,0.03089,0.000366") ? "1DOperc" :
        (damage_param == "-m b,0.024,0.36,0.1545,0.00183") ? "5DOperc" :
        (damage_param == "-m b,0.024,0.36,0.3089,0.00367") ? "10DOperc" :
        (damage_param == "-m b,0.024,0.36,0.4634,0.0055") ? "15DOperc" :
        (damage_param == "-m b,0.024,0.36,0.6179,0.00732") ? "20DOperc" :
        (damage_param == "-m b,0.024,0.36,0.772,0.0091") ? "25DOperc" :
        (damage_param == "-m b,0.024,0.36,0.926,0.011") ? "30DOperc" :  null
    when :
         readL =
        (read_length_param == "-ld Uni,20,49") ? "ReadL20-49":
        (read_length_param == "-ld Uni,50,69") ? "ReadL50-69":
        (read_length_param == "-ld Uni,70,89") ? "ReadL70-89": null

    script:

    gargammel_path = "/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/gargammel/bin/gargammel"
    data_path = "/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/aDNA_paper/dataset/Reference_genomes/Yersinia_pestis.fna"
    Bam_name = bam_file.simpleName.replaceAll('_merged\\.fastq\\.bam$', '')
    output_prefix = "Simulation_${Bam_name}"

    """
    file_name=\$(basename "$bam_file")
    #bam_file_name2="/mnt/aDNA_paper/dataset/Mappings/\${file_name}"
    real_path=\$(readlink -f $bam_file)
    bam_file_name2=\$(echo "\$real_path" | sed 's|/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN|/mnt|')
    #bam_file_name2="/mnt/aDNA_paper/dataset/Mappings/\${real_path}"
    bam_file_name=\$(basename ${bam_file} .bam)
    misincorporation_file="results_"\${bam_file_name}"/misincorporation.txt"
    if [ "${damage_level}" == "normDOperc" ]; then
        module load PDC apptainer
        echo "Running mapDamage with observed damage level ..." &>> gargammelRun.log
        singularity exec -B /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN:/mnt /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/mapdamage2_2.2.2.sif mapDamage -i \${bam_file_name2} -r /mnt/aDNA_paper/dataset/Reference_genomes/${params.species}.fna  &>> gargammelRun.log
        cp \${misincorporation_file} . &>> gargammelRun.log 
        find . -type d -name "results_Yersinia_pestis_*" -exec rm -r {} +
    fi
    echo "Running gargammel with damage level ${damage_level}..." &>> gargammelRun.log 
    #${gargammel_path} -n 1500000 --comp 0,0,1 ${damage_param} double --loc 4.106487474 --scale 0.358874723 -o "${output_prefix}_${damage_level}" ${data_path} &>> gargammelRun.log 
    #Remove useless files
    #find . -type f \\( -name "*a.fa.gz" -o -name "*b.fa.gz" -o -name "*c.fa.gz" -o -name "*e.fa.gz" -o -name "*s1.fq.gz" -o -name "*s2.fq.gz" -o -name "*.bam.bai" \\) -delete
    /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/My_conda/bin/ngsngs -t ${params.threads} -i ${data_path} -qs 40 -r 2000000 ${read_length_param} --format fq.gz -seq PE ${damage_param} -o "${output_prefix}_${readL}_${damage_level}" &>> gargammelRun.log
    echo "Running fastp to merge simulated reads"
    #fastp -i "${output_prefix}_${readL}_${damage_level}_R1.fq.gz"  -I "${output_prefix}_${readL}_${damage_level}_R2.fq.gz" -m --merged_out "${output_prefix}_${readL}_${damage_level}.fastq.gz" fastq -w ${params.threads} &>> gargammelRun.log
    SeqPrep -f "${output_prefix}_${readL}_${damage_level}_R1.fq.gz" -r "${output_prefix}_${readL}_${damage_level}_R2.fq.gz" -s "${output_prefix}_${readL}_${damage_level}.fastq.gz" -1 "${output_prefix}_${readL}_${damage_level}_out_R1.fq.gz" -2 "${output_prefix}_${readL}_${damage_level}_out_R2.fq.gz" &>> gargammelRun.log
    echo "Changing fastq to fasta"
    zcat "${output_prefix}_${readL}_${damage_level}.fastq.gz" | seqtk seq -A | pigz -p ${params.threads} > "${output_prefix}_${readL}_${damage_level}.fa.gz" 
    echo "Removing useless files"
    find . -type f \\( -name "fastp.json" -o -name "fastp.html" -o -name "*_R1.fq.gz" -o -name "*_R2.fq.gz" -o -name "*.fastq.gz" \\) -delete 
    """
   
}


process ExtractReads {
    errorStrategy 'ignore'
    
    publishDir params.output_gargammel, mode: 'symlink'

    input:
    val bam_file

    output:
    path ("*_d.fa.gz"), emit: read_file

    script:

    def output_prefix = "Simulation_${bam_file.simpleName.replaceAll('_merged\\.fastq\\.bam$', '')}"

    """
    samtools1_21="/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/My_conda/bin/samtools"
    \${samtools1_21} fasta -@ {params.threads} -F 4  ${bam_file} | gzip >  "${output_prefix}_d.fa.gz" 
    """
}



process DIVERGENCE {

    errorStrategy 'ignore'

    publishDir params.output_gargammel, mode: 'symlink'  

    conda "${params.conda_env2}"

    maxForks 499

    input:
    tuple val(assembly), val(divergence)

    output:
    path ("*.fna"), emit: assembly

    when : 
        divergence_level = 
        (divergence == 0.00) ? "div0perc" :
        //(divergence == 0.01) ? "div1perc" :
        (divergence == 0.02) ? "div2perc" :
        //(divergence == 0.03) ? "div3perc" :
        (divergence == 0.04) ? "div4perc" :
        //(divergence == 0.05) ? "div5perc" :
        //(divergence == 0.06) ? "div6perc" :
        (divergence == 0.07) ? "div7perc" :
        //(divergence == 0.08) ? "div8perc" :
        (divergence == 0.09) ? "div9perc" :
        //(divergence == 0.10) ? "div10perc" :
        //(divergence == 0.11) ? "div11perc" :
        (divergence == 0.12) ? "div12perc" :
        //(divergence == 0.13) ? "div13perc" :
        (divergence == 0.14) ? "div14perc" :
        //(divergence == 0.15) ? "div15perc" :
        (divergence == 0.16) ? "div16perc" :
        //(divergence == 0.17) ? "div17perc" :
        //(divergence == 0.18) ? "div18perc" :
        //(divergence == 0.19) ? "div19perc" :
        (divergence == 0.20) ? "div20perc" : null
        //(divergence == 0.21) ? "div21perc" :
        //(divergence == 0.22) ? "div22perc" :
        //(divergence == 0.25) ? "div25perc" : null

    script:

    //Add divergence
    """
    echo "Running divergence script"
    # Generate reads with divergence 
    base_name=\$(basename ${assembly} .fna)
    python3 /cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/TOOLS/add_divergence_to_assembly.py --input_fasta ${params.ref_genomes}${assembly}.fna --output_fasta \${base_name}_${divergence_level}.fna --divergence ${divergence} &>> AddDIvergence.log
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
    ln -s ${genome} . 
    base_name=\$(basename ${genome} .fna)
    bowtie2-build \${base_name}.fna \${base_name} &>> indexGenomes.log
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
    ln -s ${genome} . 
    base_name=\$(basename ${genome} .fna)
    bwa index \${base_name}.fna &>> indexGenomes.log
    """
}


process divergenceRunBowtie2 {

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
    def base_name_only = fq_file.simpleName.replaceAll('_d\\.fa\\.gz$', '') // Remove the suffix to get the base name
    //def ref_species = Indexed_genomes.simpleName // Extract reference species name
    //def ref_species =  basename(Indexed_genomes)
    //def ref_species = file(Indexed_genomes).simpleName
    //samtools1_21="/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/My_conda/bin/samtools"

    // Run Bowtie2
    """
    echo "Running divergence script and Bowtie2 alignment for ${fq_file}..."

    samtools1_21="/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/My_conda/bin/samtools"

    ref_species=\$(basename ${Indexed_genomes})
    base_name=\$(basename ${fq_file} .fa.gz)

    # Map the reads to the reference using Bowtie2
    bowtie2 --quiet --very-sensitive --threads ${params.threads} -x ${Indexed_genomes} -f -U ${fq_file} | \
    \${samtools1_21} view -b -F 4 -q 1 -@ ${params.threads} | \
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
    def base_name_only = fq_file.simpleName.replaceAll('_d\\.fa\\.gz$', '') // Remove the suffix to get the base name
    //samtools1_21="/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/My_conda/bin/samtools"
    //run BWA
    """
    echo "Running divergence script and BWA alignment for ${fq_file}..."

    samtools1_21="/cfs/klemming/projects/supr/snic2022-6-144/BENJAMIN/My_conda/bin/samtools"

    ref_species=\$(basename ${Indexed_genomes})
    base_name=\$(basename ${fq_file} .fa.gz)

    # Map the reads to the reference using BWA
    bwa aln -l 16500 -n 0.01 -o 2 -t ${params.threads} ${Indexed_genomes}.fna ${fq_file} > "\${base_name}_ref-\${ref_species}_ProgBWA.sai"
    bwa samse ${Indexed_genomes}.fna "\${base_name}_ref-\${ref_species}_ProgBWA.sai" ${fq_file} | \
    \${samtools1_21} view -b -q 1 -F 4 -@ ${params.threads} | \${samtools1_21} sort -@ ${params.threads} -O bam -o "\${base_name}_ref-\${ref_species}_ProgBWA.bam" &>> divergenceRun.log

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
    damage_params = Channel.of(" ", "-m b,0.024,0.36,0.03089,0.000366", "-m b,0.024,0.36,0.1545,0.00183", "-m b,0.024,0.36,0.3089,0.00367", "-m b,0.024,0.36,0.4634,0.0055", "-m b,0.024,0.36,0.6179,0.00732", "-m b,0.024,0.36,0.772,0.0091", "-m b,0.024,0.36,0.926,0.011")


    // Create a channel from BAM files 
    bam_files_non_udg = Channel.fromPath("${params.output_bam}/${params.species}*none*.bam")

    // Read length distributions 
    read_length_param = Channel.of("-ld Uni,20,49","-ld Uni,50,69","-ld Uni,70,89")

    data_gargammel = damage_params.combine(bam_files_non_udg).combine(read_length_param)

    // Run gargammelRun for non-UDG samples
    gargammelRun(data_gargammel)
  

    divergences = Channel.of(0.00,0.02,0.04,0.07,0.09,0.12,0.14,0.16,0.20)

    divergence_files = Channel.from(params.genome_list).combine(divergences)
    divergence_files.view { "divergence file: $it" }
    DIVERGENCE(divergence_files)

    // Create a channel from the genome list
    Genome_list = DIVERGENCE.out.assembly.flatten()
    
    // Index the genomes
    indexGenomesBowtie2(Genome_list)
    indexGenomesBWA(Genome_list)
    
    //gargammelRun.out.view()

    fq_file_non_udg = gargammelRun.out.gargammel_read_file.flatten()
    
    Indexed_genomes = indexGenomes.out.flatten().collect { it.toString().split('\\.')[0] }
    indexGenomesBowtie2.out.view()
    indexGenomesBWA.out.view()
    Indexed_genomesBowtie2 = indexGenomesBowtie2.out.flatten().map { it.toString().split('\\.')[0] }
    Indexed_genomesBWA = indexGenomesBWA.out.flatten().map { "${it.toString().split('\\.')[0]}" }

    sizes = Channel.of(1000000, 100000, 10000, 5000, 1000, 500, 100)

    divergence_data_non_udg_Bowtie2 = fq_file_non_udg.combine(Indexed_genomesBowtie2).combine(sizes)

    divergence_data_non_udg_BWA = fq_file_non_udg.combine(Indexed_genomesBWA).combine(sizes)

    divergence_data_BWA = divergence_data_udg_BWA.mix(divergence_data_non_udg_BWA)
    divergence_data_Bowtie2 = divergence_data_udg_Bowtie2.mix(divergence_data_non_udg_Bowtie2)
    

    // Run divergenceRun for output from ExtractReads and gargammelRun
    divergenceRunBowtie2(divergence_data_non_udg_Bowtie2)
    divergenceRunBWA(divergence_data_non_udg_BWA)
}
