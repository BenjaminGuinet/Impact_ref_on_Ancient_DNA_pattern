nextflow.enable.dsl=2

// ============================================================
// Species included in this study (reviewer comments #1, #2,
// #18, #19 — second reviewer):
//
//   Yersinia_pestis       — bacterium  (dsDNA)
//   Parvovirus_B19        — ssDNA virus (Parvoviridae)
//   Hepatitis_B_virus     — dsDNA/ssDNA virus (Hepadnaviridae)
//   Hepatitis_C_virus     — RNA virus, ssRNA(+) (Flaviviridae)
//   Torque_teno_virus     — ssDNA virus (Anelloviridae)
//
// NOTE: RNA and ssDNA viruses evolve faster than bacteria,
// giving a wider natural divergence range and making the
// divergence-correction approach (section 3) more relevant.
// ============================================================

// ============================================================
// REVISION — aligner settings (reviewer comment on l.113-135)
//
//   (a) Bowtie2 seed mismatches: the "--sensitive" preset implies
//       -N 0 (no mismatch tolerated inside the seed), whereas BWA
//       aln with -l 16500 has its seeding switched off. Bowtie2
//       is now run with "--sensitive -N 1" (params.bowtie2_params)
//       so that one mismatch is allowed in the seed. -N given
//       explicitly overrides the preset value.
//
//   (b) Mapping-quality filter: alignments were previously kept
//       at MAPQ >= 1; the threshold is now params.mapq_threshold
//       (= 25) for both aligners. 25 is the highest bwa aln
//       threshold that still keeps uniquely placed reads at the
//       -n mismatch ceiling (bwa aln gives them MAPQ 25, 37 below
//       the ceiling, 0 to multi-mappers).
//
//   (c) MAPQ sweep (exploratory). Each run now aligns once with no
//       MAPQ filter (all mapped reads >= 30 bp) and derives the main
//       BAM at params.mapq_threshold from it. For the task at
//       subsample depth params.mapq_sweep_depth (default 1,000,000),
//       one BAM per threshold in params.mapq_sweep (default
//       0,20,25,30) is also written and published to
//       params.output_mapq_sweep_bam, together with a small
//       *_mapq_sweep_counts.tsv. These BAMs are NOT used downstream.
//       Only one task per (reads, reference, aligner) writes them, so
//       the 7 subsample depths do not produce 7 identical copies.
//       Restrict to a subset of read files with
//       --mapq_sweep_pattern '<regex>', disable with --mapq_sweep none.
//
//   (d) The number of reads aligned before any MAPQ filter is added
//       to every *_mapping_stats.tsv as the last column
//       (n_aligned_mapq0).
//
//   All settings are written into every *_mapping_stats.tsv
//   (columns mapq_threshold, aligner_params) and therefore into
//   alignment_stats_summary.tsv.
// ============================================================

// MAPQ threshold applied with "samtools view -q" after alignment.
// Defined in the .config; this fallback only applies if it is missing.
params.mapq_threshold = (params.mapq_threshold != null) ? params.mapq_threshold : 25

// MAPQ sweep (see revision note (c)). Defined in the .config; fallbacks here.
params.mapq_sweep         = (params.mapq_sweep != null)         ? params.mapq_sweep         : [0, 20, 25, 30]
params.mapq_sweep_depth   = (params.mapq_sweep_depth != null)   ? params.mapq_sweep_depth   : 1_000_000
params.mapq_sweep_pattern = (params.mapq_sweep_pattern != null) ? params.mapq_sweep_pattern : null

// Fail fast on a malformed --mapq_sweep instead of inside every task
// (errorStrategy 'ignore' would otherwise drop those tasks silently).
mapqSweepList()

// ============================================================
// Helper functions
// ============================================================

def damageLabel(String param) {
    switch (param.trim()) {
        case ""                                  : return "0DOperc"
        case "-m b,0.024,0.36,0.03089,0.000366" : return "1DOperc"
        case "-m b,0.024,0.36,0.1545,0.00183"   : return "5DOperc"
        case "-m b,0.024,0.36,0.3089,0.00367"   : return "10DOperc"
        case "-m b,0.024,0.36,0.4634,0.0055"    : return "15DOperc"
        case "-m b,0.024,0.36,0.6179,0.00732"   : return "20DOperc"
        case "-m b,0.024,0.36,0.772,0.0091"     : return "25DOperc"
        case "-m b,0.024,0.36,0.926,0.011"      : return "30DOperc"
        default: return "unknownDamage"
    }
}

def readLengthLabel(String param) {
    switch (param) {
        case "-ld Uni,30,49" : return "ReadL30-49"
        case "-ld Uni,50,69" : return "ReadL50-69"
        case "-ld Uni,70,89" : return "ReadL70-89"
        default: return "unknownReadL"
    }
}

def divergenceLabel(Number div) {
    return "div${(div * 100).toInteger()}perc"
}

// MAPQ thresholds of the exploratory sweep, as a list of strings.
// Accepts a list from the config or a comma-separated string from the CLI.
def mapqSweepList() {
    def raw = params.mapq_sweep
    // A bare "--mapq_sweep" (or "--mapq_sweep ''") is parsed by Nextflow as
    // Boolean true, so any Boolean means "no sweep".
    if (raw == null || raw instanceof Boolean) return []
    def str = raw.toString().trim()
    if (str == '' || str.toLowerCase() in ['none', 'false', 'off']) return []
    def items = (raw instanceof List) ? raw.collect { it.toString().trim() }
                                      : str.split(',').collect { it.trim() }
    items = items.findAll { it }
    def bad = items.findAll { !(it ==~ /\d+/) }
    if (bad) throw new IllegalArgumentException("params.mapq_sweep must be a list of integers (or 'none'); got: ${bad}")
    return items
}

// True for the one task per (reads, reference, aligner) that writes the sweep.
def doMapqSweep(String fq_name, n_aligned_subsample) {
    if (!mapqSweepList()) return false
    if ((n_aligned_subsample as Long) != (params.mapq_sweep_depth as Long)) return false
    def pat = params.mapq_sweep_pattern
    if (pat == null || pat instanceof Boolean || pat.toString().trim() == '') return true
    return (fq_name =~ pat.toString()).find()
}

// Bash block writing one BAM (+ index) per sweep threshold from all_mapped.bam,
// plus a counts table. Expects base_name, ref_species and samtools_bin to be set.
def mapqSweepBash(String prog, String aligner, List qs) {
    return """
    # ── MAPQ sweep (exploratory; these BAMs are not used downstream) ─────────
    sweep_base="\${base_name}_ref-\${ref_species}_${prog}"
    printf "sample\\treference\\taligner\\tmapq_min\\tn_reads\\n" > "\${sweep_base}_mapq_sweep_counts.tsv"
    for q in ${qs.join(' ')}; do
        \${samtools_bin} view -b -q \$q all_mapped.bam > "\${sweep_base}_mq\${q}.bam"
        \${samtools_bin} index "\${sweep_base}_mq\${q}.bam"
        n_q=\$(\${samtools_bin} view -c "\${sweep_base}_mq\${q}.bam")
        printf "%s\\t%s\\t%s\\t%s\\t%s\\n" "\${base_name}" "\${ref_species}" "${aligner}" "\$q" "\$n_q" \\
            >> "\${sweep_base}_mapq_sweep_counts.tsv"
    done
    """
}

def bamDirForSpecies(String species) {
    def overrides = params.output_bam_per_species ?: [:]
    return overrides.containsKey(species)
        ? overrides[species]
        : "${params.output_bam}/${species}"
}


// ============================================================
// Output folder for final generated BAM files
// ============================================================
//
// The final subsampled sorted BAM files created by the Bowtie2/BWA
// divergence runs are published here as symbolic links.
// Override from the command line, for example:
//   --output_generated_bam /vol/projects/bguinet/aDNA_damage_simulations/Simulations/generated_bam_links
params.output_generated_bam = params.output_generated_bam ?: "generated_bam_links"

// Exploratory MAPQ-sweep BAMs (revision note (c)) are published here, kept
// apart from the main generated BAMs.
params.output_mapq_sweep_bam = params.output_mapq_sweep_bam ?: "${params.output_generated_bam}/mapq_sweep"

// ============================================================
// PROCESS: ngsNgsRun
//   Simulate aDNA reads → SeqPrep merge → FASTA.
//   Keeps only merged .fa.gz + sim stats TSV.
//   All intermediate FASTQ files deleted immediately.
// ============================================================
process ngsNgsRun {

    tag "${species}_${damageLabel(damage_param)}_${readLengthLabel(read_length_param)}"

    errorStrategy 'ignore'

    publishDir params.output_ngsngs, mode: 'symlink'

    conda params.conda_env_ngsngs

    input:
    tuple val(damage_param), val(species), val(read_length_param)

    output:
    path "*.fa.gz",         emit: simulated_reads
    path "*_sim_stats.tsv", emit: sim_stats

    script:
    def damage_lbl  = damageLabel(damage_param)
    def readl_lbl   = readLengthLabel(read_length_param)
    def out_prefix  = "Simulation_${species}_${readl_lbl}_${damage_lbl}"
    def damage_flag = (damage_param.trim() == "") ? "" : damage_param

    """
    # ── 1. Simulate paired-end reads ────────────────────────────────────────
    ngsngs \\
        -t ${params.threads_ngsngs} \\
        -i ${params.ref_genomes[species]} \\
        -qs 40 \\
        -r 2000000 \\
        ${read_length_param} \\
        --format fq.gz \\
        -seq PE \\
        ${damage_flag} \\
        -o "${out_prefix}" > /dev/null 2>&1

    # ── 2. Merge with SeqPrep ────────────────────────────────────────────────
    SeqPrep \\
        -f "${out_prefix}_R1.fq.gz" \\
        -r "${out_prefix}_R2.fq.gz" \\
        -s "${out_prefix}.fastq.gz" \\
        -1 "${out_prefix}_out_R1.fq.gz" \\
        -2 "${out_prefix}_out_R2.fq.gz" > /dev/null 2>&1

    # ── 3. FASTQ → FASTA ────────────────────────────────────────────────────
    zcat "${out_prefix}.fastq.gz" \\
        | seqtk seq -A \\
        | pigz -p ${params.threads_ngsngs} > "${out_prefix}.fa.gz"

    # ── 4. Sim stats TSV ─────────────────────────────────────────────────────
    n_sim_req=2000000
    n_sim_merged=\$(zcat "${out_prefix}.fa.gz" | grep -c "^>" || true)
    printf "species\tdamage_level\tread_length\tn_simulated_requested\tn_simulated_merged\n" \\
        > "${out_prefix}_sim_stats.tsv"
    printf "${species}\t${damage_lbl}\t${readl_lbl}\t\${n_sim_req}\t\${n_sim_merged}\n" \\
        >> "${out_prefix}_sim_stats.tsv"

    # ── 5. Delete all intermediates ──────────────────────────────────────────
    rm -f "${out_prefix}_R1.fq.gz" "${out_prefix}_R2.fq.gz" \\
          "${out_prefix}_out_R1.fq.gz" "${out_prefix}_out_R2.fq.gz" \\
          "${out_prefix}.fastq.gz"
    """
}


// ============================================================
// PROCESS: ExtractReads
//   Extract aligned reads from an empirical BAM into FASTA.
// ============================================================
process ExtractReads {

    tag "${bam_file.simpleName}"

    errorStrategy 'ignore'

    publishDir params.output_ngsngs, mode: 'symlink'

    input:
    tuple val(species), path(bam_file)

    output:
    path "*_d.fa.gz", emit: read_file

    script:
    def out_prefix = "${species}_empirical_${bam_file.simpleName.replaceAll('_merged\\.fastq\\.bam$', '')}"
    """
    ${params.samtools_path} fasta \\
        -@ ${params.threads} \\
        -F 4 \\
        ${bam_file} \\
        | gzip > "${out_prefix}_d.fa.gz"
    """
}


// ============================================================
// PROCESS: DIVERGENCE
//   Introduce nucleotide divergence into a reference assembly.
//   Log deleted after completion.
// ============================================================
process DIVERGENCE {

    tag "${assembly}_${divergenceLabel(divergence as Double)}"

    errorStrategy 'ignore'

    publishDir params.output_ngsngs, mode: 'symlink'

    conda params.conda_env_tools

    input:
    tuple val(assembly), val(divergence)

    output:
    path "*.fna", emit: assembly

    script:
    def div_lbl = divergenceLabel(divergence as Double)
    """
    ${params.python_bin} ${params.add_divergence} \\
        -i  ${params.ref_genomes[assembly]} \\
        -m  ${params.msa[assembly]} \\
        -o  ${assembly}_${div_lbl}.fna \\
        -d  ${divergence} > /dev/null 2>&1
    """
}


// ============================================================
// PROCESS: indexGenomesBowtie2
//   Build Bowtie2 index. No publishDir — index only lives in
//   the work dir and is referenced by path from downstream.
// ============================================================
process indexGenomesBowtie2 {

    tag "${genome.simpleName}"

    input:
    path genome

    output:
    path "*.rev.1.bt2{,l}", emit: bt2_index

    script:
    """
    bowtie2-build ${genome} ${genome.simpleName} > /dev/null 2>&1
    """
}


// ============================================================
// PROCESS: indexGenomesBWA
//   Build BWA index. No publishDir.
// ============================================================
process indexGenomesBWA {

    tag "${genome.simpleName}"

    input:
    path genome

    output:
    path "*.fna.sa", emit: bwa_index

    script:
    """
    bwa index ${genome} > /dev/null 2>&1
    """
}


// ============================================================
// PROCESS: divergenceRunBowtie2
//   Map reads → subsample → DamageProfiler → PyDamage.
//
//   Space strategy:
//     - Full-alignment BAM is still deleted immediately after subsampling
//     - Final sorted subsampled BAM and BAI are kept
//     - Final BAM and BAI are published as symbolic links in params.output_generated_bam
//     - PyDamage PNG not generated (--pl flag kept for CSV only)
//     - DamagePlot PDF not kept
//     - Logs discarded (redirected to /dev/null)
// ============================================================
process divergenceRunBowtie2 {

    tag "${fq_file.simpleName}_n${n_aligned_subsample}"

    errorStrategy 'ignore'

    publishDir params.output_ngsngs, mode: 'symlink', pattern: "*.{csv,txt,tsv}"
    publishDir params.output_generated_bam, mode: 'symlink', pattern: "*_sorted.bam*"
    publishDir params.output_mapq_sweep_bam, mode: 'symlink', pattern: "*_mq*.bam*"

    conda params.conda_env_tools

    input:
    tuple path(fq_file), val(indexed_genome), val(n_aligned_subsample)

    output:
    path "Pydamage_*.csv",                  emit: pydamage_csv
    path "*_3p_freq_misincorporations.txt", emit: freq_3p
    path "*_5p_freq_misincorporations.txt", emit: freq_5p
    path "*_mapping_stats.tsv",             emit: mapping_stats
    path "*_sorted.bam",                    emit: sorted_bam
    path "*_sorted.bam.bai",                emit: sorted_bam_bai
    path "*_ProgBowtie2_mq*.bam",           optional: true, emit: mapq_sweep_bam
    path "*_ProgBowtie2_mq*.bam.bai",       optional: true, emit: mapq_sweep_bai
    path "*_mapq_sweep_counts.tsv",         optional: true, emit: mapq_sweep_counts

    script:
    def sweep_block = doMapqSweep(fq_file.name, n_aligned_subsample)
        ? mapqSweepBash("ProgBowtie2", "Bowtie2", mapqSweepList())
        : "# MAPQ sweep not written by this task (only at subsample depth ${params.mapq_sweep_depth})"
    """
    samtools_bin="${params.samtools_path}"
    ref_species=\$(basename ${indexed_genome})
    base_name=\$(basename ${fq_file} .fa.gz)
    out_base="\${base_name}_ref-\${ref_species}_naligned${n_aligned_subsample}_ProgBowtie2"
    sorted_bam="\${out_base}_sorted.bam"

    # ── 1. Count input reads ─────────────────────────────────────────────────
    n_input=\$(zcat ${fq_file} | grep -c "^>" || true)

    # ── 2. Align ─────────────────────────────────────────────────────────────
    #   params.bowtie2_params = "--sensitive -N 1": one mismatch allowed in
    #   the seed (the preset alone implies -N 0).
    #   All mapped reads >= 30 bp are kept first, whatever their MAPQ, so that
    #   the MAPQ-0 count and the MAPQ-sweep BAMs can be derived from them.
    bowtie2 \\
        --quiet ${params.bowtie2_params} \\
        --threads ${params.threads_bowtie2} \\
        -x ${indexed_genome} \\
        -f -U ${fq_file} \\
        | \${samtools_bin} view -b -F 4 -m 30 -@ ${params.threads_bowtie2} \\
        | \${samtools_bin} sort  -@ ${params.threads_bowtie2} -O bam \\
              -o all_mapped.bam > /dev/null 2>&1

    n_aligned_mapq0=\$(\${samtools_bin} view -c all_mapped.bam)

    #   Main alignment set (subsampling, DamageProfiler, PyDamage, AMBER):
    #   MAPQ >= params.mapq_threshold. Filtering a sorted BAM keeps it sorted.
    \${samtools_bin} view -b -q ${params.mapq_threshold} -@ ${params.threads_bowtie2} all_mapped.bam \\
        > "\${base_name}_ref-\${ref_species}_full.bam"

${sweep_block}
    rm -f all_mapped.bam

    # ── 3. Mapping rate ──────────────────────────────────────────────────────
    n_aligned=\$(\${samtools_bin} view -c "\${base_name}_ref-\${ref_species}_full.bam")
    if [ "\${n_input}" -gt 0 ]; then
        mapping_rate=\$(LC_ALL=C awk 'BEGIN {printf "%.2f", '"\${n_aligned}"'/'"\${n_input}"'*100}')
    else
        mapping_rate="NA"
    fi

    # ── 4. Subsample → sorted BAM; delete full BAM immediately ──────────────
    \${samtools_bin} view -H "\${base_name}_ref-\${ref_species}_full.bam" > header.sam
    \${samtools_bin} view  "\${base_name}_ref-\${ref_species}_full.bam" \\
        | shuf -n ${n_aligned_subsample} > random_reads.sam
    cat header.sam random_reads.sam \\
        | \${samtools_bin} view -bS - \\
        | \${samtools_bin} sort  -@ ${params.threads_bowtie2} -O bam \\
              -o "\${sorted_bam}" > /dev/null 2>&1
    rm -f "\${base_name}_ref-\${ref_species}_full.bam" header.sam random_reads.sam

    \${samtools_bin} index -@ ${params.threads_bowtie2} "\${sorted_bam}" > /dev/null 2>&1

    # ── 5. Mapping stats TSV ─────────────────────────────────────────────────
    #   Original seven columns kept in place; MAPQ threshold and aligner
    #   parameters appended so the settings travel with the results.
    printf "sample\treference\taligner\tn_input_reads\tn_total_aligned\tmapping_rate_pct\tn_subsampled\tmapq_threshold\taligner_params\tn_aligned_mapq0\n" \\
        > "\${out_base}_mapping_stats.tsv"
    printf "\${base_name}\t\${ref_species}\tBowtie2\t\${n_input}\t\${n_aligned}\t\${mapping_rate}\t${n_aligned_subsample}\t${params.mapq_threshold}\t${params.bowtie2_params}\t\${n_aligned_mapq0}\n" \\
        >> "\${out_base}_mapping_stats.tsv"

    # ── 6. DamageProfiler ────────────────────────────────────────────────────
    java -jar ${params.tools_path}/DamageProfiler-1.1-java11.jar \\
        -i "\${sorted_bam}" \\
        -o DamageProfiler_output > /dev/null 2>&1 || true
    if [ -f DamageProfiler_output/5p_freq_misincorporations.txt ]; then
        mv DamageProfiler_output/5p_freq_misincorporations.txt "\${out_base}_5p_freq_misincorporations.txt"
        mv DamageProfiler_output/3p_freq_misincorporations.txt "\${out_base}_3p_freq_misincorporations.txt"
    else
        printf "pos\tA\tC\tG\tT\n" > "\${out_base}_5p_freq_misincorporations.txt"
        printf "pos\tA\tC\tG\tT\n" > "\${out_base}_3p_freq_misincorporations.txt"
    fi
    rm -rf DamageProfiler_output

    # ── 7. PyDamage (CSV only, no PNG) ──────────────────────────────────────
    pydamage analyze "\${sorted_bam}" -p ${params.threads_bowtie2} -g > /dev/null 2>&1 || true
    if [ -f pydamage_results/pydamage_results.csv ]; then
        cp pydamage_results/pydamage_results.csv "Pydamage_\${out_base}.csv"
    else
        printf "reference,coverage,qvalue,predicted_accuracy,null_model_p4c,damage_model_p4c\n" \\
            > "Pydamage_\${out_base}.csv"
    fi
    rm -rf pydamage_results

    # ── 8. Keep final BAM + BAI ───────────────────────────────────────────────
    # Do not delete "\${sorted_bam}" or "\${sorted_bam}.bai";
    # Nextflow publishes them as symlinks to params.output_generated_bam.
    """
}


// ============================================================
// PROCESS: divergenceRunBWA
//   Identical logic to divergenceRunBowtie2 but using BWA
//   aln/samse, appropriate for short ancient DNA fragments.
// ============================================================
process divergenceRunBWA {

    tag "${fq_file.simpleName}_n${n_aligned_subsample}"

    errorStrategy 'ignore'

    publishDir params.output_ngsngs, mode: 'symlink', pattern: "*.{csv,txt,tsv}"
    publishDir params.output_generated_bam, mode: 'symlink', pattern: "*_sorted.bam*"
    publishDir params.output_mapq_sweep_bam, mode: 'symlink', pattern: "*_mq*.bam*"

    conda params.conda_env_tools

    input:
    tuple path(fq_file), val(indexed_genome), val(n_aligned_subsample)

    output:
    path "Pydamage_*.csv",                  emit: pydamage_csv
    path "*_3p_freq_misincorporations.txt", emit: freq_3p
    path "*_5p_freq_misincorporations.txt", emit: freq_5p
    path "*_mapping_stats.tsv",             emit: mapping_stats
    path "*_sorted.bam",                    emit: sorted_bam
    path "*_sorted.bam.bai",                emit: sorted_bam_bai
    path "*_ProgBWA_mq*.bam",               optional: true, emit: mapq_sweep_bam
    path "*_ProgBWA_mq*.bam.bai",           optional: true, emit: mapq_sweep_bai
    path "*_mapq_sweep_counts.tsv",         optional: true, emit: mapq_sweep_counts

    script:
    def sweep_block = doMapqSweep(fq_file.name, n_aligned_subsample)
        ? mapqSweepBash("ProgBWA", "BWA", mapqSweepList())
        : "# MAPQ sweep not written by this task (only at subsample depth ${params.mapq_sweep_depth})"
    """
    samtools_bin="${params.samtools_path}"
    ref_species=\$(basename ${indexed_genome})
    base_name=\$(basename ${fq_file} .fa.gz)
    out_base_bwa="\${base_name}_ref-\${ref_species}_naligned${n_aligned_subsample}_ProgBWA"
    sorted_bam="\${out_base_bwa}_sorted.bam"

    # ── 1. Count input reads ─────────────────────────────────────────────────
    n_input=\$(zcat ${fq_file} | grep -c "^>" || true)

    # ── 2. Align with BWA aln/samse ──────────────────────────────────────────
    bwa aln ${params.bwa_aln_params} -t ${params.threads_bwa} \\
        ${indexed_genome}.fna ${fq_file} \\
        > "\${base_name}_ref-\${ref_species}.sai" 2> /dev/null

    #   All mapped reads >= 30 bp are kept first, whatever their MAPQ, so that
    #   the MAPQ-0 count and the MAPQ-sweep BAMs can be derived from them.
    bwa samse ${indexed_genome}.fna \\
        "\${base_name}_ref-\${ref_species}.sai" ${fq_file} \\
        | \${samtools_bin} view -b -F 4 -m 30 -@ ${params.threads_bwa} \\
        | \${samtools_bin} sort  -@ ${params.threads_bwa} -O bam \\
              -o all_mapped.bam > /dev/null 2>&1
    rm -f "\${base_name}_ref-\${ref_species}.sai"

    n_aligned_mapq0=\$(\${samtools_bin} view -c all_mapped.bam)

    #   Main alignment set (subsampling, DamageProfiler, PyDamage, AMBER):
    #   MAPQ >= params.mapq_threshold. Filtering a sorted BAM keeps it sorted.
    \${samtools_bin} view -b -q ${params.mapq_threshold} -@ ${params.threads_bwa} all_mapped.bam \\
        > "\${base_name}_ref-\${ref_species}_full.bam"

${sweep_block}
    rm -f all_mapped.bam

    # ── 3. Mapping rate ──────────────────────────────────────────────────────
    n_aligned=\$(\${samtools_bin} view -c "\${base_name}_ref-\${ref_species}_full.bam")
    if [ "\${n_input}" -gt 0 ]; then
        mapping_rate=\$(LC_ALL=C awk 'BEGIN {printf "%.2f", '"\${n_aligned}"'/'"\${n_input}"'*100}')
    else
        mapping_rate="NA"
    fi

    # ── 4. Subsample → sorted BAM; delete full BAM immediately ──────────────
    \${samtools_bin} view -H "\${base_name}_ref-\${ref_species}_full.bam" > header.sam
    \${samtools_bin} view  "\${base_name}_ref-\${ref_species}_full.bam" \\
        | shuf -n ${n_aligned_subsample} > random_reads.sam
    cat header.sam random_reads.sam \\
        | \${samtools_bin} view -bS - \\
        | \${samtools_bin} sort  -@ ${params.threads_bwa} -O bam \\
              -o "\${sorted_bam}" > /dev/null 2>&1
    rm -f "\${base_name}_ref-\${ref_species}_full.bam" header.sam random_reads.sam

    \${samtools_bin} index -@ ${params.threads_bwa} "\${sorted_bam}" > /dev/null 2>&1

    # ── 5. Mapping stats TSV ─────────────────────────────────────────────────
    printf "sample\treference\taligner\tn_input_reads\tn_total_aligned\tmapping_rate_pct\tn_subsampled\tmapq_threshold\taligner_params\tn_aligned_mapq0\n" \\
        > "\${out_base_bwa}_mapping_stats.tsv"
    printf "\${base_name}\t\${ref_species}\tBWA\t\${n_input}\t\${n_aligned}\t\${mapping_rate}\t${n_aligned_subsample}\t${params.mapq_threshold}\t${params.bwa_aln_params}\t\${n_aligned_mapq0}\n" \\
        >> "\${out_base_bwa}_mapping_stats.tsv"

    # ── 6. DamageProfiler ────────────────────────────────────────────────────
    java -jar ${params.tools_path}/DamageProfiler-1.1-java11.jar \\
        -i "\${sorted_bam}" \\
        -o DamageProfiler_output > /dev/null 2>&1 || true
    if [ -f DamageProfiler_output/5p_freq_misincorporations.txt ]; then
        mv DamageProfiler_output/5p_freq_misincorporations.txt "\${out_base_bwa}_5p_freq_misincorporations.txt"
        mv DamageProfiler_output/3p_freq_misincorporations.txt "\${out_base_bwa}_3p_freq_misincorporations.txt"
    else
        printf "pos\tA\tC\tG\tT\n" > "\${out_base_bwa}_5p_freq_misincorporations.txt"
        printf "pos\tA\tC\tG\tT\n" > "\${out_base_bwa}_3p_freq_misincorporations.txt"
    fi
    rm -rf DamageProfiler_output

    # ── 7. PyDamage (CSV only, no PNG) ──────────────────────────────────────
    pydamage analyze "\${sorted_bam}" -p ${params.threads_bwa} -g > /dev/null 2>&1 || true
    if [ -f pydamage_results/pydamage_results.csv ]; then
        cp pydamage_results/pydamage_results.csv "Pydamage_\${out_base_bwa}.csv"
    else
        printf "reference,coverage,qvalue,predicted_accuracy,null_model_p4c,damage_model_p4c\n" \\
            > "Pydamage_\${out_base_bwa}.csv"
    fi
    rm -rf pydamage_results

    # ── 8. Keep final BAM + BAI ───────────────────────────────────────────────
    # Do not delete "\${sorted_bam}" or "\${sorted_bam}.bai";
    # Nextflow publishes them as symlinks to params.output_generated_bam.
    """
}



process computeReferenceBiasAmber {

    tag "${out_base}"

    publishDir params.output_reference_bias, mode: 'copy', pattern: "*.{tsv,pdf}"

    conda params.conda_env_tools

    input:
    tuple val(out_base), path(sorted_bam), path(sorted_bai)

    output:
    path "${out_base}_mismatch_by_read_length.tsv", emit: mismatch_by_read_length
    path "${out_base}_reference_bias.summary.tsv",  emit: refbias_summary
    path "${out_base}_amber.pdf", optional: true,   emit: amber_pdf

    script:
    """
    set -euo pipefail

    amber_prefix="${out_base}_amber"

    # Make sure pysam finds the BAM index locally
    if [ ! -e "${sorted_bam}.bai" ]; then
        ln -sf "${sorted_bai}" "${sorted_bam}.bai" || true
    fi

    ${params.python_bin} ${params.amber_single} \\
        --bam ${sorted_bam} \\
        --sample ${out_base} \\
        --output "\${amber_prefix}"

    # Normalize mismatch-by-read-length output.
    # Prefer native AMBER table if your amber_single.py writes it.
    if [ -s "\${amber_prefix}.mismatch_by_read_length.tsv" ]; then
        cp "\${amber_prefix}.mismatch_by_read_length.tsv" "${out_base}_mismatch_by_read_length.tsv"
    else
        ${params.python_bin} - <<'PY'
import csv

amber_prefix = "${out_base}_amber"
out_base = "${out_base}"
bam_path = "${sorted_bam}"

inp = amber_prefix + ".tables.tsv"
out = "${out_base}_mismatch_by_read_length.tsv"

with open(out, "w", newline="") as fo:
    w = csv.writer(fo, delimiter="\\t")
    w.writerow([
        "sample",
        "bam_path",
        "read_length",
        "n_reads",
        "mean_mismatches_per_bp",
        "ci95_mean_mismatches_per_bp",
        "pooled_mismatches_per_bp",
    ])

    with open(inp, newline="") as fi:
        r = csv.DictReader(fi, delimiter="\\t")
        for row in r:
            if row.get("table") != "mismatch_by_read_length":
                continue

            ci95 = "NA"
            extra = row.get("extra", "")
            if extra.startswith("ci95="):
                ci95 = extra.replace("ci95=", "")

            w.writerow([
                out_base,
                bam_path,
                row.get("x", "NA"),
                row.get("n", "0"),
                row.get("value", "NA"),
                ci95,
                row.get("value", "NA"),
            ])
PY
    fi

    ${params.python_bin} - <<'PY'
import csv, math

inp = "${out_base}_mismatch_by_read_length.tsv"
out = "${out_base}_reference_bias.summary.tsv"

rows = []
with open(inp, newline="") as f:
    r = csv.DictReader(f, delimiter="\\t")
    for row in r:
        try:
            rl = float(row["read_length"])
            n = float(row["n_reads"])
            mean = float(row["mean_mismatches_per_bp"])
            pooled = float(row["pooled_mismatches_per_bp"])
        except Exception:
            continue
        if n > 0:
            rows.append((rl, n, mean, pooled))

header = [
    "sample",
    "refbias_n_read_lengths",
    "refbias_total_reads_by_length",
    "refbias_min_read_length",
    "refbias_max_read_length",
    "refbias_weighted_mean_mismatches_per_bp",
    "refbias_weighted_pooled_mismatches_per_bp",
    "refbias_slope_mismatch_per_bp_per_read_bp",
    "refbias_delta_predicted_max_minus_min",
]

def fmt(x):
    if x is None or math.isnan(x) or math.isinf(x):
        return "NA"
    return f"{x:.10g}"

if rows:
    total_n = sum(n for rl, n, mean, pooled in rows)
    min_rl = min(rl for rl, n, mean, pooled in rows)
    max_rl = max(rl for rl, n, mean, pooled in rows)
    weighted_mean = sum(mean * n for rl, n, mean, pooled in rows) / total_n
    weighted_pooled = sum(pooled * n for rl, n, mean, pooled in rows) / total_n

    xbar = sum(rl * n for rl, n, mean, pooled in rows) / total_n
    ybar = weighted_mean
    ssx = sum(n * (rl - xbar) ** 2 for rl, n, mean, pooled in rows)

    if ssx > 0:
        slope = sum(n * (rl - xbar) * (mean - ybar) for rl, n, mean, pooled in rows) / ssx
        delta = slope * (max_rl - min_rl)
    else:
        slope = float("nan")
        delta = float("nan")

    values = [
        "${out_base}",
        str(len(rows)),
        str(int(total_n)),
        fmt(min_rl),
        fmt(max_rl),
        fmt(weighted_mean),
        fmt(weighted_pooled),
        fmt(slope),
        fmt(delta),
    ]
else:
    values = [
        "${out_base}",
        "0",
        "0",
        "NA",
        "NA",
        "NA",
        "NA",
        "NA",
        "NA",
    ]

with open(out, "w", newline="") as f:
    w = csv.writer(f, delimiter="\\t")
    w.writerow(header)
    w.writerow(values)
PY

    if [ -s "\${amber_prefix}.pdf" ]; then
        cp "\${amber_prefix}.pdf" "${out_base}_amber.pdf"
    fi
    """
}

// ============================================================
// PROCESS: computeDamageMetrics
//   Compute stochasticity, intensity, PyDamage summary, and
//   per-position mismatch tables into one TSV + damage plot.
//   Only the TSV and PDF remain in the work dir after completion.
// ============================================================
process computeDamageMetrics {

    tag "${out_base}"

    errorStrategy 'ignore'

    publishDir params.output_damage_metrics, mode: 'copy'

    conda params.conda_env_tools

    input:
    tuple val(out_base),
          path(freq_5p_sample),
          path(freq_3p_sample),
          path(freq_5p_ref),
          path(freq_3p_ref),
          path(pydamage_csv),
          path(stats_tsv)

    output:
    path "${out_base}_damage_metrics.tsv", emit: metrics_tsv
    path "${out_base}_damage_plot.pdf",    emit: metrics_pdf

    script:
    """
    ${params.python_bin} ${params.compute_damage_metrics} \\
        --sample_5p ${freq_5p_sample} \\
        --sample_3p ${freq_3p_sample} \\
        --ref_5p    ${freq_5p_ref} \\
        --ref_3p    ${freq_3p_ref} \\
        --pydamage  ${pydamage_csv} \\
        --stats     ${stats_tsv} \\
        --out       ${out_base}_damage_metrics.tsv \\
        --plot      ${out_base}_damage_plot.pdf
    """
}


// ============================================================
// PROCESS: collectDamageMetrics
//   Merge all per-job damage_metrics TSVs into one summary.
//   Files staged into metrics/ subdir to avoid name collisions.
// ============================================================
process collectDamageMetrics {

    publishDir params.output_damage_metrics, mode: 'copy'

    input:
    path 'metrics/*'

    output:
    path "damage_metrics_summary.tsv"

    script:
    """
    head -n 1 \$(ls metrics/*_damage_metrics.tsv | head -1) \\
        > damage_metrics_summary.tsv
    for f in metrics/*_damage_metrics.tsv; do
        tail -n +2 "\$f" >> damage_metrics_summary.tsv
    done
    """
}


// ============================================================
// PROCESS: collectStats
//   Merge all per-job mapping_stats TSVs into one summary.
// ============================================================
process collectStats {

    publishDir params.output_ngsngs, mode: 'copy'

    input:
    path stats_files

    output:
    path "alignment_stats_summary.tsv"

    script:
    """
    printf "sample\treference\taligner\tn_input_reads\tn_total_aligned\tmapping_rate_pct\tn_subsampled\tmapq_threshold\taligner_params\tn_aligned_mapq0\n" \\
        > alignment_stats_summary.tsv
    for f in ${stats_files}; do
        tail -n +2 "\$f" >> alignment_stats_summary.tsv
    done
    """
}


// ============================================================
// WORKFLOW
// ============================================================
workflow {

    species_list = Channel.from(params.genome_list)

    // ── Damage parameter grid ────────────────────────────────────────────────
    damage_params = Channel.of(
        " ",                                 // 0 %  — no damage
        "-m b,0.024,0.36,0.03089,0.000366",  //  1 %
        "-m b,0.024,0.36,0.1545,0.00183",    //  5 %
        "-m b,0.024,0.36,0.3089,0.00367",    // 10 %
        "-m b,0.024,0.36,0.4634,0.0055",     // 15 %
        "-m b,0.024,0.36,0.6179,0.00732",    // 20 %
        "-m b,0.024,0.36,0.772,0.0091",      // 25 %
        "-m b,0.024,0.36,0.926,0.011"        // 30 %
    )

    // ── Read-length distribution grid ────────────────────────────────────────
    read_length_params = Channel.of(
      "-ld Uni,30,49",
      "-ld Uni,50,69",
      "-ld Uni,70,89"
    )

    sim_input = damage_params
        .combine(species_list)
        .combine(read_length_params)

    ngsNgsRun(sim_input)

    // ── Empirical reads ──────────────────────────────────────────────────────
    bam_files = species_list.flatMap { sp ->
        def dir = bamDirForSpecies(sp)
        file("${dir}/${sp}*none*.bam").collect { bam -> tuple(sp, bam) }
    }
    ExtractReads(bam_files)

    external_reads = params.input_fastq
        ? Channel.fromPath(params.input_fastq)
        : Channel.empty()

    all_reads = ngsNgsRun.out.simulated_reads
        .flatten()
        .mix(ExtractReads.out.read_file.flatten())
        .mix(external_reads)

    // ── Divergence grid ──────────────────────────────────────────────────────
    divergences = Channel.of(0.00, 0.02, 0.04, 0.07, 0.09, 0.12, 0.14, 0.16, 0.20)
    //divergences = Channel.of(0.00, 0.12)

    divergence_inputs = species_list.combine(divergences)
    DIVERGENCE(divergence_inputs)

    // Tag each diverged assembly with its species for same-species pairing
    diverged_assemblies = DIVERGENCE.out.assembly.flatten()
        .map { fna ->
            def sp = fna.simpleName.replaceAll(/_div\d+perc$/, '')
            tuple(sp, fna)
        }

    indexGenomesBowtie2(diverged_assemblies.map { sp, fna -> fna })
    indexGenomesBWA(diverged_assemblies.map { sp, fna -> fna })

    bt2_prefixes = indexGenomesBowtie2.out.bt2_index
        .flatten()
        .map { f ->
            def prefix = f.toString().replaceAll(/\.rev\.1\.bt2l?$/, '')
            def sp     = new File(prefix).name.replaceAll(/_div\d+perc$/, '')
            tuple(sp, prefix)
        }

    bwa_prefixes = indexGenomesBWA.out.bwa_index
        .flatten()
        .map { f ->
            def prefix = f.toString().replaceAll(/\.fna\.sa$/, '')
            def sp     = new File(prefix).name.replaceAll(/_div\d+perc$/, '')
            tuple(sp, prefix)
        }

    // Tag each read file with its species
    // Simulated reads: filename starts with Simulation_<species>_
    // Empirical reads: filename starts with <species>_empirical_
    tagged_reads = all_reads.map { fq ->
        def sp = params.genome_list.find { s ->
            fq.name.startsWith("Simulation_${s}_") || fq.name.startsWith("${s}_empirical_")
        }
        tuple(sp, fq)
    }

    // ── Subsampling depths ───────────────────────────────────────────────────
    n_aligned_subsample = Channel.of(1_000_000, 100_000, 10_000, 5_000, 1_000, 500, 100)
    //n_aligned_subsample = Channel.of(1_000_000, 5_000)

    // combine(by: 0) joins on species — reads only map to same-species references
    data_bowtie2 = tagged_reads
        .combine(bt2_prefixes, by: 0)
        .combine(n_aligned_subsample)
        .map { sp, fq, prefix, n -> tuple(fq, prefix, n) }

    data_bwa = tagged_reads
        .combine(bwa_prefixes, by: 0)
        .combine(n_aligned_subsample)
        .map { sp, fq, prefix, n -> tuple(fq, prefix, n) }

    divergenceRunBowtie2(data_bowtie2)
    divergenceRunBWA(data_bwa)

    // ── Collect mapping stats ────────────────────────────────────────────────
    all_mapping_stats = divergenceRunBowtie2.out.mapping_stats
        .mix(divergenceRunBWA.out.mapping_stats)
        .collect()

    collectStats(all_mapping_stats)

    // ── Damage metrics ───────────────────────────────────────────────────────
    // Join key: strip _naligned{N}_Prog{X}_{5p|3p}_freq_misincorporations.txt
    def make_key = { f ->
        f.name.replaceAll(/_naligned\d+_Prog\w+_(5p|3p)_freq_misincorporations\.txt$/, '')
    }
    def get_n = { f ->
        def m = (f.name =~ /_naligned(\d+)_/)[0]
        m ? m[1].toLong() : 0L
    }

    def metrics_channel = { freq_5p_ch, freq_3p_ch, pydamage_ch, stats_ch ->

        t5     = freq_5p_ch.flatten().map { f -> tuple(make_key(f), get_n(f), f) }
        t3     = freq_3p_ch.flatten().map { f -> tuple(make_key(f), get_n(f), f) }
        paired = t5.join(t3, by: [0, 1])

        tpyd = pydamage_ch.flatten().map { f ->
            def base = f.name.replaceAll(/^Pydamage_/, '')
            def key  = base.replaceAll(/_naligned\d+_Prog\w+\.csv$/, '')
            def m    = (base =~ /_naligned(\d+)_/)[0]
            def n    = m ? m[1].toLong() : 0L
            tuple(key, n, f)
        }

        paired_pyd = paired.join(tpyd, by: [0, 1])

        ref_1M = paired_pyd
            .filter { it[1] == 1_000_000L }
            .map    { key, n, f5, f3, fpyd -> tuple(key, f5, f3) }

        subs = paired_pyd
            .filter { it[1] != 1_000_000L }
            .map    { key, n, f5, f3, fpyd ->
                def out_base = f5.name.replace('_5p_freq_misincorporations.txt', '')
                tuple(key, out_base, f5, f3, fpyd)
            }

        // Stats are matched on (sample, subsampling depth) so that each
        // computeDamageMetrics job receives the stats of its own subsample.
        stats_tagged = stats_ch.flatten().map { f ->
            def key = f.name.replaceAll(/_naligned\d+_Prog\w+_mapping_stats\.tsv$/, '')
            tuple(key, get_n(f), f)
        }

        subs
            .map  { key, out_base, f5, f3, fpyd ->
                tuple(key, get_n(f5), out_base, f5, f3, fpyd)
            }
            .join(stats_tagged, by: [0, 1])
            .map  { key, n, out_base, f5, f3, fpyd, stats ->
                tuple(key, out_base, f5, f3, fpyd, stats)
            }
            .combine(ref_1M, by: 0)
            .map { key, out_base, f5_sub, f3_sub, fpyd_sub, stats, f5_ref, f3_ref ->
                tuple(out_base, f5_sub, f3_sub, f5_ref, f3_ref, fpyd_sub, stats)
            }
    }

    metrics_bt2 = metrics_channel(
        divergenceRunBowtie2.out.freq_5p,
        divergenceRunBowtie2.out.freq_3p,
        divergenceRunBowtie2.out.pydamage_csv,
        divergenceRunBowtie2.out.mapping_stats
    )
    metrics_bwa = metrics_channel(
        divergenceRunBWA.out.freq_5p,
        divergenceRunBWA.out.freq_3p,
        divergenceRunBWA.out.pydamage_csv,
        divergenceRunBWA.out.mapping_stats
    )

    def bam_key = { f ->
        def name = f.name
        name = name.replaceFirst(/\.bai$/, '')
        name = name.replaceFirst(/_sorted\.bam$/, '')
        return name
    }

    def amber_input_channel = { bam_ch, bai_ch ->
        bams = bam_ch.flatten().map { f -> tuple(bam_key(f), f) }
        bais = bai_ch.flatten().map { f -> tuple(bam_key(f), f) }
        bams.join(bais, by: 0).map { key, bam, bai -> tuple(key, bam, bai) }
    }

    amber_bt2 = amber_input_channel(
        divergenceRunBowtie2.out.sorted_bam,
        divergenceRunBowtie2.out.sorted_bam_bai
    )

    amber_bwa = amber_input_channel(
        divergenceRunBWA.out.sorted_bam,
        divergenceRunBWA.out.sorted_bam_bai
    )

    computeReferenceBiasAmber(amber_bt2.mix(amber_bwa))

    computeDamageMetrics(metrics_bt2.mix(metrics_bwa))

    collectDamageMetrics(
        computeDamageMetrics.out.metrics_tsv
            .unique { it.name }
            .collect()
    )
}
