#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * ============================================================================
 *  Empirical aDNA Damage Profiling Pipeline
 * ============================================================================
 *  Purpose of this revised version:
 *   1. Download/merge/trim empirical reads.
 *   2. Map all trimmed reads first to the closest/target reference genome.
 *   3. Extract read IDs that map to this target reference after the same filters.
 *   4. Build a target-confirmed FASTQ file.
 *   5. Remap this SAME target-confirmed read set to all references in the panel.
 *
 *  This addresses the reviewer concern that mapping raw enriched/metagenomic
 *  libraries independently to each divergent reference may generate damage
 *  profiles from different biological read pools.
 *
 *  ---------------------------------------------------------------------------
 *  CHANGES IN THIS VERSION (aligner-consistent read selection, MAPQ sweep)
 *  ---------------------------------------------------------------------------
 *  (A) Default minimum MAPQ is now 25 (params.min_mapq), the compromise chosen
 *      after comparing thresholds 1/20/25/30 on the ancient HSV SG1 sample.
 *
 *  (B) The target-confirmed read set is now defined PER ALIGNER. Previously the
 *      read set was selected with bwa aln and then remapped with both bwa aln
 *      and Bowtie2, so the Bowtie2 branch inherited bwa's sensitivity profile.
 *      Now TARGET_BWA_MAP feeds BWA_MAP and TARGET_BOWTIE2_MAP feeds BOWTIE2_MAP,
 *      each through FILTER_FASTQ_TO_TARGET_READS (which carries the aligner name
 *      in its tuple and writes one target-confirmed FASTQ per aligner).
 *
 *  (C) BWA_MAP and BOWTIE2_MAP align once, then write one BAM per threshold in
 *      params.publish_mapqs (default [1, 20, 25, 30]) next to the main BAM in
 *      the mapping results folder, named *_mq<N>.bam. Only the main BAM
 *      (filtered at params.min_mapq) continues downstream; the others are for
 *      exploratory use.
 *
 *  (D) Bowtie2 runs with --sensitive -N 1 (one mismatch allowed in the seed),
 *      which is the usual setting for short damaged reads.
 *
 *  DOWNLOAD_FASTQ, MERGE_FASTQ, FASTP, BOWTIE2_INDEX and BWA_INDEX are
 *  unchanged, so they stay cached under -resume. Everything from the target
 *  mapping onwards will rerun.
 *
 *  ---------------------------------------------------------------------------
 *  EARLIER CHANGES (read-name / mate-pair bug fix)
 *  ---------------------------------------------------------------------------
 *  (1) MERGE_FASTQ
 *      Old behaviour detected mates with `ls *_1.fastq.gz`, so files named
 *      *_R1.fastq.gz fell through to the "OTHER" branch and were concatenated
 *      into a single R1. FASTP then took its single-end path and both mates of
 *      every fragment survived as independent reads. Unpaired files were also
 *      appended to R1, which breaks fastp's lock-step pairing.
 *      New behaviour pairs files by stripping the mate token (_1/_2 or _R1/_R2),
 *      verifies the mate exists, and emits paired files and genuinely single-end
 *      files separately (_merged_R1 / _merged_R2 / _merged_SE).
 *
 *  (2) FASTP
 *      Runs the PE collapse on the paired files and a separate SE trim on the
 *      unpaired files, then concatenates both. Critically, it now normalises
 *      read names: only the first whitespace-delimited token is kept and any
 *      trailing /1 or /2 is removed. bwa strips /1 and /2 when writing QNAME,
 *      so without this the BAM names and the FASTQ IDs can never match.
 *
 *  (3) FILTER_FASTQ_TO_TARGET_READS
 *      Two hard guards replace the old "fail only at exactly zero" check:
 *        - duplicate QNAMEs in the target BAM (mates not collapsed upstream)
 *        - fewer than 99% of target read IDs recovered from the FASTQ
 *      Either condition now stops the pipeline instead of silently producing a
 *      tiny read set that propagates all the way to the metrics table.
 *  ---------------------------------------------------------------------------
 *
 *  Optional: your sample sheet can contain a column called target_ref_name.
 *  If target_ref_name is absent, the workflow infers the target reference from
 *  singlegenome_species. For example:
 *
 *  sample_name  singlegenome_species  archive_accession
 *  Ypestis      Yersinia_pestis       ERR2862147
 *
 *  will automatically select a reference whose basename is either exactly
 *  Yersinia_pestis or starts with Yersinia_pestis__, e.g.
 *  Yersinia_pestis__GCF_000009065.1.fna.
 *
 *  You can still override this by adding target_ref_name explicitly, e.g.
 *  target_ref_name = Yersinia_pestis__GCF_000009065.1
 * ============================================================================
 */

// Optional local FASTQ reuse. Set in nextflow.config, for example:
// params.local_fastq_dir  = '/vol/projects/bguinet/aDNA_damage_simulations/Empiric/Fastqs'
// params.local_fastq_only = true   // optional: fail if local FASTQs are missing instead of downloading
params.local_fastq_dir  = params.local_fastq_dir ?: ""
params.local_fastq_only = params.local_fastq_only ?: false

// Minimum MAPQ applied at every mapping step (TARGET_BWA_MAP, TARGET_BOWTIE2_MAP,
// BWA_MAP, BOWTIE2_MAP). The main BAM that feeds subsampling, DamageProfiler
// and the metrics uses this threshold. Override with --min_mapq N.
params.min_mapq = params.min_mapq ?: 25

// Additional MAPQ thresholds for which BWA_MAP and BOWTIE2_MAP also write a
// BAM (named *_mq<N>.bam) in the mapping results folder. These BAMs are NOT
// used downstream; they exist for exploratory comparisons. Override with
// --publish_mapqs 1,20,25,30 (comma-separated).
params.publish_mapqs = params.publish_mapqs ?: [1, 20, 25, 30]

// ── Processes ───────────────────────────────────────────────────────────────

process DOWNLOAD_FASTQ {
    tag "${sample_name}_${accession}"
    publishDir "${params.outdir}/${sample_name}/fastq", mode: 'symlink'

    input:
    tuple val(species), val(sample_name), val(target_ref_name), val(accession)

    output:
    tuple val(species), val(sample_name), val(target_ref_name), path("*.fastq.gz"), emit: fastq

    script:
    """
    #!/bin/bash
    set -euo pipefail

    ACC="${accession}"
    SAMPLE="${sample_name}"
    TARGET="${target_ref_name}"
    LOCAL_DIR="${params.local_fastq_dir}"
    LOCAL_ONLY="${params.local_fastq_only}"

    echo "[INFO] Sample:      \$SAMPLE"
    echo "[INFO] Accession:   \$ACC"
    echo "[INFO] Target ref:  \$TARGET"

    # ── Optional local FASTQ reuse ──────────────────────────────────────────
    # If params.local_fastq_dir is set, first look for files containing the
    # accession (preferred), then sample name, then target reference prefix.
    # Example matches for ACC=ERR2862147:
    #   Yersinia_pestis_ERR2862147_1.fastq.gz
    #   Yersinia_pestis_ERR2862147_2.fastq.gz
    if [[ -n "\$LOCAL_DIR" ]]; then
        if [[ ! -d "\$LOCAL_DIR" ]]; then
            echo "[ERROR] params.local_fastq_dir does not exist: \$LOCAL_DIR" >&2
            exit 1
        fi

        echo "[INFO] Looking for local FASTQs in: \$LOCAL_DIR"
        shopt -s nullglob

        LOCAL_FILES=( "\$LOCAL_DIR"/*"\$ACC"*.fastq.gz "\$LOCAL_DIR"/*"\$ACC"*.fq.gz )

        if [[ \${#LOCAL_FILES[@]} -eq 0 ]]; then
            LOCAL_FILES=( "\$LOCAL_DIR"/*"\$SAMPLE"*.fastq.gz "\$LOCAL_DIR"/*"\$SAMPLE"*.fq.gz )
        fi

        if [[ \${#LOCAL_FILES[@]} -eq 0 ]]; then
            LOCAL_FILES=( "\$LOCAL_DIR"/*"\$TARGET"*.fastq.gz "\$LOCAL_DIR"/*"\$TARGET"*.fq.gz )
        fi

        if [[ \${#LOCAL_FILES[@]} -gt 0 ]]; then
            echo "[INFO] Found \${#LOCAL_FILES[@]} local FASTQ file(s); linking instead of downloading:"
            for f in "\${LOCAL_FILES[@]}"; do
                echo "  - \$f"
                b=\$(basename "\$f")
                # Normalise .fq.gz to .fastq.gz so downstream MERGE_FASTQ sees it.
                if [[ "\$b" == *.fq.gz ]]; then
                    b="\${b%.fq.gz}.fastq.gz"
                fi
                ln -sf "\$f" "\$b"
            done
        else
            echo "[WARN] No local FASTQs found for accession/sample/target: \$ACC / \$SAMPLE / \$TARGET"
            if [[ "\$LOCAL_ONLY" == "true" || "\$LOCAL_ONLY" == "1" || "\$LOCAL_ONLY" == "yes" ]]; then
                echo "[ERROR] local_fastq_only=true, so refusing to download missing FASTQs." >&2
                exit 1
            fi
        fi
    fi

    # If local FASTQs were linked, skip downloading. Otherwise download normally.
    if ls *.fastq.gz 1>/dev/null 2>&1; then
        echo "[INFO] Using local FASTQ file(s):"
        ls -lh *.fastq.gz
    else
        echo "[INFO] No local FASTQs available; downloading accession \$ACC"

        # ── Case 1: Full URL (http/https/ftp) — download directly with wget ──
        if [[ "\$ACC" == http://* ]] || [[ "\$ACC" == https://* ]] || [[ "\$ACC" == ftp://* ]]; then
            FNAME=\$(basename "\$ACC")
            echo "[INFO] Full URL detected — downloading \$FNAME directly with wget"
            wget -q "\$ACC" -O "\$FNAME"

            # Gzip if not already compressed
            if [[ "\$FNAME" != *.gz ]]; then
                pigz -p ${task.cpus} "\$FNAME"
            fi

        # ── Case 2: ERS/ERP/PRJEB/SAMEA — resolve via ENA API ──
        elif [[ "\$ACC" == ERS* ]] || [[ "\$ACC" == ERP* ]] || [[ "\$ACC" == PRJEB* ]] || [[ "\$ACC" == SAMEA* ]]; then
            echo "[INFO] Resolving sample accession \$ACC via ENA API..."
            curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=\${ACC}&result=read_run&fields=run_accession,fastq_ftp&format=tsv" \
                > ena_report.tsv
            cat ena_report.tsv

            tail -n +2 ena_report.tsv | while IFS=\$'\t' read -r run_acc ftp_urls; do
                if [[ -z "\$ftp_urls" ]]; then
                    echo "[WARN] No FTP URLs for \$run_acc, trying fasterq-dump..."
                    fasterq-dump --split-files --threads ${task.cpus} --temp . "\$run_acc"
                    for f in \${run_acc}*.fastq; do pigz -p ${task.cpus} "\$f"; done
                else
                    IFS=';' read -ra URLS <<< "\$ftp_urls"
                    for url in "\${URLS[@]}"; do
                        url=\$(echo "\$url" | sed 's/^ *//;s/ *\$//')
                        fname=\$(basename "\$url")
                        echo "[INFO] Downloading \$fname"
                        wget -q "ftp://\$url" -O "\$fname" || wget -q "http://\$url" -O "\$fname"
                    done
                fi
            done

        # ── Case 3: Standard run accession (ERR/SRR/DRR/CRR etc.) — fasterq-dump ──
        else
            echo "[INFO] Run accession \$ACC — using fasterq-dump"
            fasterq-dump --split-files --threads ${task.cpus} --temp . "\$ACC"
            for f in \${ACC}*.fastq; do pigz -p ${task.cpus} "\$f"; done
        fi
    fi

    # Verify output
    if ! ls *.fastq.gz 1>/dev/null 2>&1; then
        echo "[ERROR] No fastq.gz files for \$ACC" >&2
        exit 1
    fi

    echo "[INFO] Final FASTQ outputs:"
    ls -lh *.fastq.gz
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// FIXED: pair detection.
//
// Mates are identified by stripping a trailing _1/_2 or _R1/_R2 token and
// checking that the partner file actually exists. Files with no partner are
// kept apart as genuine single-end input instead of being appended to R1,
// which used to break fastp's lock-step pairing.
//
// Outputs (any subset may be produced):
//   ${sample_name}_merged_R1.fastq.gz / _merged_R2.fastq.gz   paired
//   ${sample_name}_merged_SE.fastq.gz                         unpaired
// ─────────────────────────────────────────────────────────────────────────────
process MERGE_FASTQ {
    tag "${sample_name}"

    input:
    tuple val(species), val(sample_name), val(target_ref_name), path(fastqs)

    output:
    tuple val(species), val(sample_name), val(target_ref_name), path("${sample_name}_merged_*.fastq.gz"), emit: merged

    script:
    """
    #!/bin/bash
    set -euo pipefail
    shopt -s nullglob

    if command -v pigz >/dev/null 2>&1; then
        COMPRESS="pigz -p ${task.cpus}"
    else
        COMPRESS="gzip -c"
    fi

    R1=(); R2=(); SE=()

    # The mate token is _1/_2 or _R1/_R2, optionally followed by a lane or
    # chunk suffix such as _001 (bcl2fastq style: SAMP_S1_L001_R1_001.fastq.gz).
    for f in *.fastq.gz *.fq.gz; do

        [[ "\$f" == ${sample_name}_merged_* ]] && continue
        [[ "\$f" == split_*_[12].fastq.gz ]]   && continue
        [[ "\$f" == orphan_*.fastq.gz ]]       && continue

        if [[ "\$f" =~ ^(.*)_(R?)1(_[0-9]+)?\\.(fastq|fq)\\.gz\$ ]]; then
            stem="\${BASH_REMATCH[1]}"; rtag="\${BASH_REMATCH[2]}"
            tail="\${BASH_REMATCH[3]}"; ext="\${BASH_REMATCH[4]}"
            mate="\${stem}_\${rtag}2\${tail}.\${ext}.gz"
            if [[ -s "\$mate" ]]; then
                R1+=( "\$f" ); R2+=( "\$mate" )
            else
                echo "[WARN] No mate found for \$f — treating as single-end"
                SE+=( "\$f" )
            fi

        elif [[ "\$f" =~ ^(.*)_(R?)2(_[0-9]+)?\\.(fastq|fq)\\.gz\$ ]]; then
            stem="\${BASH_REMATCH[1]}"; rtag="\${BASH_REMATCH[2]}"
            tail="\${BASH_REMATCH[3]}"; ext="\${BASH_REMATCH[4]}"
            mate="\${stem}_\${rtag}1\${tail}.\${ext}.gz"
            if [[ ! -s "\$mate" ]]; then
                echo "[WARN] Orphan R2 \$f — treating as single-end"
                SE+=( "\$f" )
            fi
            # paired R2 files are picked up through their R1 above

        else
            SE+=( "\$f" )
        fi
    done

    # Some archives ship BOTH mates inside ONE file, told apart only by the
    # /1 and /2 suffix on the read name. This is common for ".mapped" extracts
    # pulled back out of a BAM, where each mate was kept independently: the two
    # mate counts are NOT equal and the file is not cleanly interleaved, so the
    # pairing has to be done by read name rather than by position.
    #
    # Left alone, fastp would trim these as independent single-end reads and
    # both mates of every fragment would enter the analysis, inflating every
    # count ~2x and making the subsample depths meaningless.
    KEEP=()
    nsplit=0
    if [[ \${#SE[@]} -gt 0 ]]; then
        for f in "\${SE[@]}"; do

            NMATE=\$(zcat "\$f" | awk 'NR%4==1 { split(\$0, a, " "); id=a[1]; t=substr(id, length(id)-1); if (t=="/1" || t=="/2") c++ } END { print c+0 }')
            if [[ "\$NMATE" -eq 0 ]]; then
                KEEP+=( "\$f" )
                continue
            fi

            nsplit=\$(( nsplit + 1 ))
            echo "[INFO] \$f carries /1 and /2 mate suffixes — pairing by read name"

            # One record per line: base_id, mate, then the 4 FASTQ lines.
            # Sorting by (base_id, mate) is what guarantees R1 and R2 come out
            # in the same order, which fastp requires.
            zcat "\$f" | paste - - - - | awk -F'\\t' 'BEGIN { OFS="\\t" } {
                    split(\$1, a, " "); id=a[1]; sub(/^@/, "", id);
                    t=substr(id, length(id)-1);
                    if (t=="/1" || t=="/2") { base=substr(id, 1, length(id)-2); mate=substr(t, 2, 1) }
                    else                    { base=id; mate="0" }
                    print base, mate, \$1, \$2, \$3, \$4
                }' | LC_ALL=C sort -k1,1 -k2,2 > recs_\${nsplit}.tsv

            # base ids seen with two different mate values = complete pairs
            cut -f1,2 recs_\${nsplit}.tsv | uniq | cut -f1 | uniq -d > pairids_\${nsplit}.txt

            NTOT=\$(wc -l < recs_\${nsplit}.tsv)
            NPAIR=\$(wc -l < pairids_\${nsplit}.txt)
            echo "[INFO]   \$NTOT reads -> \$NPAIR complete pairs, \$(( NTOT - 2 * NPAIR )) orphan(s)"

            TAB=\$(printf '\\t')

            if [[ "\$NPAIR" -gt 0 ]]; then
                LC_ALL=C join -t "\$TAB" -1 1 -2 1 pairids_\${nsplit}.txt recs_\${nsplit}.tsv > paired_\${nsplit}.tsv
                awk -F'\\t' '\$2=="1" { print \$3"\\n"\$4"\\n"\$5"\\n"\$6 }' paired_\${nsplit}.tsv | \$COMPRESS > split_\${nsplit}_1.fastq.gz
                awk -F'\\t' '\$2=="2" { print \$3"\\n"\$4"\\n"\$5"\\n"\$6 }' paired_\${nsplit}.tsv | \$COMPRESS > split_\${nsplit}_2.fastq.gz
                R1+=( "split_\${nsplit}_1.fastq.gz" )
                R2+=( "split_\${nsplit}_2.fastq.gz" )
                rm -f paired_\${nsplit}.tsv
            fi

            # mates with no partner stay single-end rather than being discarded
            LC_ALL=C join -t "\$TAB" -v 2 -1 1 -2 1 pairids_\${nsplit}.txt recs_\${nsplit}.tsv \\
                | awk -F'\\t' '{ print \$3"\\n"\$4"\\n"\$5"\\n"\$6 }' \\
                | \$COMPRESS > orphan_\${nsplit}.fastq.gz

            NORPH=\$(zcat orphan_\${nsplit}.fastq.gz | awk 'END { print NR/4 }')
            if [[ "\${NORPH%.*}" -gt 0 ]]; then
                KEEP+=( "orphan_\${nsplit}.fastq.gz" )
            fi

            rm -f recs_\${nsplit}.tsv
        done
    fi
    SE=( \${KEEP[@]+"\${KEEP[@]}"} )

    if [[ \${#R1[@]} -eq 0 && \${#SE[@]} -eq 0 ]]; then
        echo "[ERROR] No fastq.gz files found to merge" >&2
        exit 1
    fi

    if [[ \${#R1[@]} -gt 0 ]]; then
        echo "[INFO] Paired-end: \${#R1[@]} file pair(s)"
        for i in "\${!R1[@]}"; do
            echo "  - \${R1[\$i]}  +  \${R2[\$i]}"
        done
        cat "\${R1[@]}" > ${sample_name}_merged_R1.fastq.gz
        cat "\${R2[@]}" > ${sample_name}_merged_R2.fastq.gz
    fi

    if [[ \${#SE[@]} -gt 0 ]]; then
        echo "[INFO] Single-end: \${#SE[@]} file(s)"
        for f in "\${SE[@]}"; do echo "  - \$f"; done
        cat "\${SE[@]}" > ${sample_name}_merged_SE.fastq.gz
    fi

    echo "[INFO] Merged output files:"
    ls -lh ${sample_name}_merged_*.fastq.gz
    """
}

process FASTP {
    tag "${sample_name}"
    publishDir "${params.outdir}/${sample_name}/trimmed", mode: 'symlink'

    input:
    tuple val(species), val(sample_name), val(target_ref_name), path(reads)

    output:
    tuple val(species), val(sample_name), val(target_ref_name), path("${sample_name}_trimmed.fastq.gz"), emit: trimmed
    path("${sample_name}_fastp*.json"), emit: report

    script:
    def files = (reads instanceof List) ? reads : [reads]
    def r1    = files.find { it.name.endsWith('_merged_R1.fastq.gz') }
    def r2    = files.find { it.name.endsWith('_merged_R2.fastq.gz') }
    def se    = files.find { it.name.endsWith('_merged_SE.fastq.gz') }

    // An R1 with no R2 should never reach here, but if it does, trim it as SE.
    def se_in = se ?: ((r1 && !r2) ? r1 : null)

    def pe_block = (r1 && r2) ? """
    echo "[INFO] Paired-end input: collapsing overlapping mates"
    fastp -i ${r1} -I ${r2} \\
        --merge --merged_out pe_merged.fastq.gz \\
        --overlap_len_require 11 \\
        --length_required 30 \\
        --thread ${task.cpus} \\
        --json ${sample_name}_fastp_pe.json \\
        --html /dev/null \\
        2> fastp_pe.log
    PARTS+=( pe_merged.fastq.gz )
""" : '    echo "[INFO] No paired-end files for this sample"'

    def se_block = (se_in) ? """
    echo "[INFO] Single-end input: trimming"
    fastp -i ${se_in} \\
        --length_required 30 \\
        --thread ${task.cpus} \\
        --json ${sample_name}_fastp_se.json \\
        --html /dev/null \\
        -o se_trimmed.fastq.gz \\
        2> fastp_se.log
    PARTS+=( se_trimmed.fastq.gz )
""" : '    echo "[INFO] No single-end files for this sample"'

    """
    #!/bin/bash
    set -euo pipefail

    if command -v pigz >/dev/null 2>&1; then
        COMPRESS="pigz -p ${task.cpus}"
    else
        COMPRESS="gzip -c"
    fi

    PARTS=()

${pe_block}

${se_block}

    if [[ \${#PARTS[@]} -eq 0 ]]; then
        echo "[ERROR] fastp produced no output for ${sample_name}" >&2
        exit 1
    fi

    # Concatenate and normalise read names in a single pass.
    zcat "\${PARTS[@]}" \\
        | awk 'NR%4==1 { split(\$0, a, " "); id=a[1]; sub(/\\/[12]\$/, "", id); print id; next }
               { print }
               END { printf("[INFO] %d reads written\\n", NR/4) > "/dev/stderr" }' \\
        | \$COMPRESS \\
        > ${sample_name}_trimmed.fastq.gz

    echo "[INFO] First read name after normalisation:"
    zcat ${sample_name}_trimmed.fastq.gz | head -1 || true

    ls -lh ${sample_name}_trimmed.fastq.gz
    """
}

process BOWTIE2_INDEX {
    tag "${ref_name}"
    publishDir "${params.outdir}/bt2_indices/${species}", mode: 'copy'

    input:
    tuple val(species), val(ref_name), path(ref_fasta)

    output:
    tuple val(species), val(ref_name), path(ref_fasta), path("${ref_name}*.bt2"), emit: index

    script:
    """
    bowtie2-build --threads ${task.cpus} ${ref_fasta} ${ref_name}
    """
}

process BWA_INDEX {
    tag "${ref_name}"
    publishDir "${params.outdir}/bwa_indices/${species}", mode: 'copy'

    input:
    tuple val(species), val(ref_name), path(ref_fasta)

    output:
    tuple val(species), val(ref_name), path(ref_fasta), path("*.{amb,ann,bwt,pac,sa}"), emit: index

    script:
    """
    bwa index ${ref_fasta}
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// STEP 1: map all trimmed reads to the closest/target reference only.
// This defines the target-confirmed read set, once per aligner, so that the
// reads later remapped with bwa aln were selected with bwa aln, and the reads
// remapped with Bowtie2 were selected with Bowtie2.
// ─────────────────────────────────────────────────────────────────────────────
process TARGET_BWA_MAP {
    tag "${sample_name}_target_${target_ref_name}_BWA"
    publishDir "${params.outdir}/${sample_name}/target_read_filter/BWA/${target_ref_name}", mode: 'symlink'

    input:
    tuple val(species), val(sample_name), val(target_ref_name), path(reads), val(ref_name), path(ref_fasta), path(bwa_idx)

    output:
    tuple val(species), val(sample_name), val(target_ref_name), val("BWA"), path(reads), path("${sample_name}__target__${target_ref_name}__BWA.bam"), emit: target_bam

    script:
    """
    #!/bin/bash
    set -euo pipefail

    bwa aln -l 16500 -n 0.01 -o 2 -t ${task.cpus} ${ref_fasta} ${reads} > target.sai

    bwa samse ${ref_fasta} target.sai ${reads} \
    | samtools view -bS -F 4 -q ${params.min_mapq} -m 30 - \
    | samtools sort -@ ${task.cpus} -o ${sample_name}__target__${target_ref_name}__BWA.bam

    samtools index ${sample_name}__target__${target_ref_name}__BWA.bam
    """
}

process TARGET_BOWTIE2_MAP {
    tag "${sample_name}_target_${target_ref_name}_Bowtie2"
    publishDir "${params.outdir}/${sample_name}/target_read_filter/Bowtie2/${target_ref_name}", mode: 'symlink'

    input:
    tuple val(species), val(sample_name), val(target_ref_name), path(reads), val(ref_name), path(ref_fasta), path(bt2_idx)

    output:
    tuple val(species), val(sample_name), val(target_ref_name), val("Bowtie2"), path(reads), path("${sample_name}__target__${target_ref_name}__Bowtie2.bam"), emit: target_bam

    script:
    """
    #!/bin/bash
    set -euo pipefail

    bowtie2 --sensitive -N 1 -x ${ref_name} \
        -U ${reads} \
        --threads ${task.cpus} \
        --no-unal \
        2> ${sample_name}__target__${target_ref_name}__bt2.log \
    | samtools view -bS -F 4 -q ${params.min_mapq} -m 30 - \
    | samtools sort -@ ${task.cpus} -o ${sample_name}__target__${target_ref_name}__Bowtie2.bam

    samtools index ${sample_name}__target__${target_ref_name}__Bowtie2.bam
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// NEW STEP 2: extract read IDs from the target BAM and keep only those reads
// in the FASTQ. The resulting FASTQ is the same starting read set that will be
// remapped to every divergent reference.
// Requires seqkit in your environment.
//
// FIXED: two hard guards. The old check only fired when extraction returned
// exactly zero reads, which let a 312/57117 recovery propagate silently into
// the damage profiles and the metrics table.
// ─────────────────────────────────────────────────────────────────────────────
process FILTER_FASTQ_TO_TARGET_READS {
    tag "${sample_name}_target_confirmed_reads_${aligner}"
    publishDir "${params.outdir}/${sample_name}/target_read_filter/filtered_fastq/${aligner}", mode: 'symlink'

    input:
    tuple val(species), val(sample_name), val(target_ref_name), val(aligner), path(reads), path(target_bam)

    output:
    tuple val(species), val(sample_name), val(target_ref_name), val(aligner), path("${sample_name}_target_confirmed_${aligner}.fastq.gz"), emit: target_reads
    path("${sample_name}_target_read_ids_${aligner}.txt"), emit: read_ids

    script:
    def ids_file   = "${sample_name}_target_read_ids_${aligner}.txt"
    def fastq_out  = "${sample_name}_target_confirmed_${aligner}.fastq.gz"
    """
    #!/bin/bash
    set -euo pipefail

    if command -v pigz >/dev/null 2>&1; then
        COMPRESS="pigz -p ${task.cpus}"
    else
        COMPRESS="gzip -c"
    fi

    NBAM=\$(samtools view -c -F 4 -q 1 ${target_bam})

    samtools view -F 4 -q 1 ${target_bam} \\
        | awk '{print \$1}' \\
        | LC_ALL=C sort -u \\
        > ${ids_file}

    NIDS=\$(wc -l < ${ids_file})
    echo "[INFO] ${sample_name} [${aligner}]: \$NBAM target alignments, \$NIDS distinct read names (ref ${target_ref_name})"

    if [[ "\$NIDS" -eq 0 ]]; then
        echo "[ERROR] No target-confirmed reads found for ${sample_name}" >&2
        exit 1
    fi

    # ── Guard 1: duplicated QNAMEs ──────────────────────────────────────────
    # bwa samse and bowtie2 (without -k) emit exactly one record per read and
    # no secondary alignments, so NBAM should equal NIDS. A ~2:1 ratio means R1
    # and R2 of the same fragment are both in the library and collapse onto one
    # name once the aligner strips /1 and /2 — i.e. the mates were never merged
    # upstream. Every count and subsample depth downstream would be inflated
    # about twofold.
    DUP=\$(( NBAM - NIDS ))
    MAXDUP=\$(awk -v n="\$NBAM" 'BEGIN{printf "%d", n*0.01}')
    if [[ "\$DUP" -gt "\$MAXDUP" ]]; then
        echo "[ERROR] \$DUP of \$NBAM target alignments share a read name with another alignment." >&2
        echo "[ERROR] Mates were not collapsed upstream. Check MERGE_FASTQ pairing and the" >&2
        echo "[ERROR] fastp --merge step for ${sample_name}." >&2
        exit 1
    fi

    # Match FASTQ records by sequence ID, not by the full FASTQ header/name.
    # FASTP normalises names so that the BAM QNAME and the FASTQ ID are identical.
    seqkit grep -f ${ids_file} ${reads} \\
        | \$COMPRESS \\
        > ${fastq_out}

    NFASTQ=\$(zcat ${fastq_out} | awk 'END{printf "%d", NR/4}')
    echo "[INFO] ${sample_name} [${aligner}]: recovered \$NFASTQ of \$NIDS target-confirmed reads"

    # ── Guard 2: incomplete recovery ────────────────────────────────────────
    MIN=\$(awk -v n="\$NIDS" 'BEGIN{printf "%d", n*0.99}')
    if [[ "\$NFASTQ" -lt "\$MIN" ]]; then
        echo "[ERROR] Only \$NFASTQ of \$NIDS target read IDs were found in the FASTQ." >&2
        echo "[ERROR] This is a BAM/FASTQ read-name mismatch, not a biological result." >&2
        echo "[ERROR] BAM read IDs:" >&2
        head -3 ${ids_file} >&2
        echo "[ERROR] FASTQ header:" >&2
        zcat ${reads} | head -1 >&2 || true
        exit 1
    fi

    echo "[INFO] Target-confirmed FASTQ (${aligner}):"
    ls -lh ${fastq_out}
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// STEP 3: remap the target-confirmed reads to every reference in the panel.
//
// Each process aligns once, keeps every mapped read of >= 30 bp in a sorted
// "all_mapped.bam", then derives:
//   - the main BAM   (<sample>__<ref>__<aligner>.bam, MAPQ >= params.min_mapq),
//     which is the only one passed downstream;
//   - one BAM per threshold in params.publish_mapqs
//     (<sample>__<ref>__<aligner>_mq<N>.bam), published alongside for
//     exploratory comparisons and otherwise unused.
// Filtering a coordinate-sorted BAM with `samtools view -q` preserves the
// sort order, so no re-sort is needed per threshold.
// ─────────────────────────────────────────────────────────────────────────────
process BOWTIE2_MAP {
    tag "${sample_name}_vs_${ref_name}"
    publishDir "${params.outdir}/${sample_name}/mapping/Bowtie2/${ref_name}", mode: 'symlink'

    input:
    tuple val(species), val(sample_name), path(reads), val(ref_name), path(ref_fasta), path(bt2_idx)

    output:
    tuple val(species), val(sample_name), val(ref_name), val("Bowtie2"), path("${sample_name}__${ref_name}__Bowtie2.bam"), emit: bam
    path("${sample_name}__${ref_name}__Bowtie2_mq*.bam{,.bai}"), optional: true, emit: mapq_bams
    path("${sample_name}__${ref_name}__bt2.log"), emit: log

    script:
    def prefix = "${sample_name}__${ref_name}__Bowtie2"
    def publish_mapqs = (params.publish_mapqs instanceof List)
        ? params.publish_mapqs.collect { it.toString().trim() }
        : params.publish_mapqs.toString().split(',').collect { it.trim() }.findAll { it }
    """
    #!/bin/bash
    set -euo pipefail

    bowtie2 --sensitive -N 1 -x ${ref_name} \\
        -U ${reads} \\
        --threads ${task.cpus} \\
        --no-unal \\
        2> ${sample_name}__${ref_name}__bt2.log \\
    | samtools view -bS -F 4 -m 30 - \\
    | samtools sort -@ ${task.cpus} -o all_mapped.bam

    # main BAM used downstream
    samtools view -b -q ${params.min_mapq} all_mapped.bam > ${prefix}.bam
    samtools index ${prefix}.bam

    # one extra BAM per published threshold (exploratory only)
    for q in ${publish_mapqs.join(' ')}; do
        samtools view -b -q \$q all_mapped.bam > ${prefix}_mq\${q}.bam
        samtools index ${prefix}_mq\${q}.bam
        echo "[INFO] ${prefix}_mq\${q}.bam: \$(samtools view -c ${prefix}_mq\${q}.bam) reads"
    done
    echo "[INFO] ${prefix}.bam (MAPQ>=${params.min_mapq}, used downstream): \$(samtools view -c ${prefix}.bam) reads"

    rm -f all_mapped.bam
    """
}

process BWA_MAP {
    tag "${sample_name}_vs_${ref_name}"
    publishDir "${params.outdir}/${sample_name}/mapping/BWA/${ref_name}", mode: 'symlink'

    input:
    tuple val(species), val(sample_name), path(reads), val(ref_name), path(ref_fasta), path(bwa_idx)

    output:
    tuple val(species), val(sample_name), val(ref_name), val("BWA"), path("${sample_name}__${ref_name}__BWA.bam"), emit: bam
    path("${sample_name}__${ref_name}__BWA_mq*.bam{,.bai}"), optional: true, emit: mapq_bams

    script:
    def prefix = "${sample_name}__${ref_name}__BWA"
    def publish_mapqs = (params.publish_mapqs instanceof List)
        ? params.publish_mapqs.collect { it.toString().trim() }
        : params.publish_mapqs.toString().split(',').collect { it.trim() }.findAll { it }
    """
    #!/bin/bash
    set -euo pipefail

    bwa aln -l 16500 -n 0.01 -o 2 -t ${task.cpus} ${ref_fasta} ${reads} > aln.sai

    bwa samse ${ref_fasta} aln.sai ${reads} \\
    | samtools view -bS -F 4 -m 30 - \\
    | samtools sort -@ ${task.cpus} -o all_mapped.bam

    # main BAM used downstream
    samtools view -b -q ${params.min_mapq} all_mapped.bam > ${prefix}.bam
    samtools index ${prefix}.bam

    # one extra BAM per published threshold (exploratory only)
    for q in ${publish_mapqs.join(' ')}; do
        samtools view -b -q \$q all_mapped.bam > ${prefix}_mq\${q}.bam
        samtools index ${prefix}_mq\${q}.bam
        echo "[INFO] ${prefix}_mq\${q}.bam: \$(samtools view -c ${prefix}_mq\${q}.bam) reads"
    done
    echo "[INFO] ${prefix}.bam (MAPQ>=${params.min_mapq}, used downstream): \$(samtools view -c ${prefix}.bam) reads"

    rm -f all_mapped.bam
    """
}

process MAPPING_STATS {
    tag "${sample_name}_${ref_name}_${aligner}"
    publishDir "${params.outdir}/${sample_name}/mapping_stats/${aligner}/${ref_name}", mode: 'symlink'

    input:
    tuple val(species), val(sample_name), val(ref_name), val(aligner), path(bam)

    output:
    tuple val(species), val(sample_name), val(ref_name), val(aligner), path(bam), path("mapping_stats.tsv"), emit: stats

    script:
    """
    #!/vol/projects/bguinet/TOOLS/My_conda/bin/python3
    import pysam
    import numpy as np

    bam = pysam.AlignmentFile("${bam}", "rb")
    lengths, nms, n = [], [], 0

    for r in bam.fetch(until_eof=True):
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        n += 1
        lengths.append(r.query_alignment_length)
        nm = r.get_tag("NM") if r.has_tag("NM") else 0
        alen = r.query_alignment_length
        if alen > 0:
            nms.append(1.0 - nm / alen)

    bam.close()

    mean_len = np.mean(lengths) if lengths else 0
    mean_id  = np.mean(nms) * 100 if nms else 0

    with open("mapping_stats.tsv", "w") as f:
        print("sample\tspecies\tref_name\taligner\tn_total_aligned\tmean_read_length\tAvrg_Perc_ID", file=f)
        print(f"${sample_name}\t${species}\t${ref_name}\t${aligner}\t{n}\t{mean_len:.1f}\t{mean_id:.2f}", file=f)
    """.stripIndent()
}
process SUBSAMPLE_BAM {
    tag "${sample_name}_${ref_name}_${aligner}_n${n_reads}"
    publishDir "${params.outdir}/${sample_name}/subsampled/${aligner}/${ref_name}/nreads_${n_reads}", mode: 'symlink'

    input:
    tuple val(species), val(sample_name), val(ref_name), val(aligner), path(bam), path(stats), val(n_reads)

    output:
    tuple val(species), val(sample_name), val(ref_name), val(aligner), val(n_reads), path("subsampled.bam"), path(stats), emit: sub_bam

    script:
    """
    #!/bin/bash
    set -euo pipefail

    TOTAL=\$(samtools view -c ${bam})

    if [[ \$TOTAL -le ${n_reads} ]]; then
        cp ${bam} subsampled.bam
    else
        FRAC=\$(awk -v n=${n_reads} -v total=\$TOTAL 'BEGIN {printf "%.12f", n/total}')
        FRAC=\${FRAC#0}
        samtools view -@ ${task.cpus} -bs 42\${FRAC} ${bam} > subsampled.bam

        NEW_TOTAL=\$(samtools view -c subsampled.bam)
        if [[ \$NEW_TOTAL -gt ${n_reads} ]]; then
            samtools view -H subsampled.bam > subsampled.tmp.sam
            set +o pipefail
            samtools view subsampled.bam | \
            awk 'BEGIN{srand(42)} {print rand()"\t"\$0}' | \
            sort -k1,1n | \
            cut -f2- | \
            head -n ${n_reads} >> subsampled.tmp.sam
            set -o pipefail
            samtools view -b subsampled.tmp.sam | samtools sort -@ ${task.cpus} -o subsampled.tmp.bam
            rm -f subsampled.tmp.sam
            mv subsampled.tmp.bam subsampled.bam
        fi
    fi

    samtools index subsampled.bam
    """
}

process PASSTHROUGH_FULL_BAM {
    tag "${sample_name}_${ref_name}_${aligner}_all"

    input:
    tuple val(species), val(sample_name), val(ref_name), val(aligner), path(bam), path(stats)

    output:
    tuple val(species), val(sample_name), val(ref_name), val(aligner), val("all"), path("full.bam"), path(stats), emit: full_bam

    script:
    """
    #!/bin/bash
    set -euo pipefail

    cp ${bam} full.bam
    samtools index full.bam
    """
}

process DAMAGE_PROFILER {
    tag "${sample_name}_${ref_name}_${aligner}_n${n_reads}"
    publishDir "${params.outdir}/${sample_name}/damage_profiles/${aligner}/${ref_name}/nreads_${n_reads}", mode: 'symlink'

    input:
    tuple val(species), val(sample_name), val(ref_name), val(aligner), val(n_reads), path(bam), path(stats)

    output:
    tuple val(species), val(sample_name), val(ref_name), val(aligner), val(n_reads),
          path("DamageProfiler_out/5p_freq_misincorporations.txt"),
          path("DamageProfiler_out/3p_freq_misincorporations.txt"),
          path(stats),
          emit: profiles

    script:
    """
    #!/bin/bash
    set -euo pipefail

    # Re-index the BAM — the .bai is not staged from upstream processes
    samtools index ${bam}

    # Run DamageProfiler with error capture
    java -jar /vol/projects/bguinet/aDNA_damage_simulations/TOOLS/DamageProfiler-1.1-java11.jar \
        -i ${bam} \
        -o DamageProfiler_out \
        -l 100 \
        > dp_stdout.log 2> dp_stderr.log || {
            echo "[WARN] DamageProfiler exited non-zero. stderr:" >&2
            cat dp_stderr.log >&2
        }

    # Check if output was actually produced with data, not just header
    if [[ ! -f DamageProfiler_out/5p_freq_misincorporations.txt ]] || [[ \$(wc -l < DamageProfiler_out/5p_freq_misincorporations.txt) -le 1 ]]; then
        echo "[WARN] DamageProfiler produced no data rows — falling back to empty header" >&2
        if [[ -f dp_stderr.log ]]; then cat dp_stderr.log >&2; fi
        mkdir -p DamageProfiler_out
        echo -e "Pos\tC>T\tG>A\tA>C\tA>G\tA>T\tC>A\tC>G\tG>C\tG>T\tT>A\tT>C\tT>G" > DamageProfiler_out/5p_freq_misincorporations.txt
        echo -e "Pos\tC>T\tG>A\tA>C\tA>G\tA>T\tC>A\tC>G\tG>C\tG>T\tT>A\tT>C\tT>G" > DamageProfiler_out/3p_freq_misincorporations.txt
    fi
    """
}

process COMPUTE_METRICS {
    tag "${sample_name}_${ref_name}_${aligner}_n${n_reads}"
    publishDir "${params.outdir}/${sample_name}/metrics/${aligner}/${ref_name}/nreads_${n_reads}", mode: 'symlink'

    input:
    tuple val(species), val(sample_name), val(ref_name), val(aligner), val(n_reads), path(dp_5p), path(dp_3p), path(stats)
    path(scripts_dir)

    output:
    tuple val(species), path("metrics_${sample_name}__${ref_name}__${aligner}__n${n_reads}.tsv"), emit: metrics
    path("damage_plot_*.pdf"), optional: true, emit: plot

    script:
    """
    #!/bin/bash
    set +e

    export PATH=/home/bguinet/miniforge3/envs/My_conda/bin:\$PATH

    python3 ${scripts_dir}/compute_damage_metrics.py \
        --sample_5p ${dp_5p} \
        --sample_3p ${dp_3p} \
        --stats     ${stats} \
        --label     "${sample_name}__${ref_name}__${aligner}__n${n_reads}" \
        --out       metrics_${sample_name}__${ref_name}__${aligner}__n${n_reads}.tsv \
        --plot      damage_plot_${sample_name}__${ref_name}__${aligner}__n${n_reads}.pdf \
        2>&1 || true
    """
}

process MERGE_METRICS {
    tag "${species}"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(species), path(metric_files)

    output:
    path("${species}_all_metrics.tsv"), optional: true, emit: merged

    script:
    """
    #!/home/bguinet/miniforge3/envs/My_conda/bin/python3
    import pandas as pd
    import glob
    import sys

    files = glob.glob("metrics_*.tsv")

    if not files:
        pd.DataFrame().to_csv("${species}_all_metrics.tsv", sep="\t", index=False)
        sys.exit(0)

    dfs = []
    for f in sorted(files):
        try:
            dfs.append(pd.read_csv(f, sep="\t"))
        except Exception as e:
            print(f"WARN: skipping {f}: {e}", file=sys.stderr)

    if dfs:
        merged = pd.concat(dfs, ignore_index=True)

        # Parse the label: {sample_name}__{ref_name}__{aligner}__n{depth}
        # ref_name may contain __, so use regex anchored on known aligner names.
        if "sample" in merged.columns:
            pat = r'^(.+?)__(.+)__(BWA|Bowtie2)__n(.+)\$'
            extracted = merged["sample"].str.extract(pat)
            if extracted.shape[1] == 4:
                merged["sample_name"]  = extracted[0]
                merged["ref_name"]     = extracted[1]
                merged["aligner"]      = extracted[2]
                merged["n_subsampled"] = extracted[3]

        merged.to_csv("${species}_all_metrics.tsv", sep="\t", index=False)
    else:
        pd.DataFrame().to_csv("${species}_all_metrics.tsv", sep="\t", index=False)
    """.stripIndent()
}

// ── Workflow ────────────────────────────────────────────────────────────────

workflow {

    // ── Parse sample sheet ──────────────────────────────────────────────────
    // Required columns:
    // sample_name, singlegenome_species, archive_accession
    // Optional column:
    // target_ref_name, target_reference, or closest_reference
    ch_samples = Channel
        .fromPath(params.sample_sheet)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def sp = (row.singlegenome_species ?: "").trim()

            if (!sp) {
                error "Could not read column 'singlegenome_species'. Your sample sheet must be TAB-separated and contain at least: sample_name, singlegenome_species, archive_accession. Check with: cat -A ${params.sample_sheet}"
            }

            def species_key = ""

            if (sp.contains("Yersinia")) species_key = "Yersinia"
            else if (sp.contains("Hepatitis")) species_key = "Hepatitis"
            else if (sp.contains("Mycobacterium")) species_key = "Mycobacterium"
            else if (sp.contains("Variola")) species_key = "Variola"
            else if (sp.contains("herpesvirus") || sp.contains("Herpesvirus") || sp.contains("alphaherpesvirus")) species_key = "Herpesvirus"
            else if (sp.contains("Salmonella")) species_key = "Salmonella"
            else if (sp.contains("Parvovirus") || sp.contains("parvovirus")) species_key = "Parvovirus"
            else species_key = sp.replaceAll(/\s+/, "_")

            // Target-reference selection:
            // If target_ref_name is provided in the sample sheet, use it.
            // Otherwise, infer it from singlegenome_species. For example:
            //   singlegenome_species = Yersinia_pestis
            // will match a reference basename such as:
            //   Yersinia_pestis__GCF_000009065.1
            def target_ref_name = (row.target_ref_name ?: row.target_reference ?: row.closest_reference ?: "").trim()
            if (!target_ref_name) {
                target_ref_name = sp.replaceAll(/\s+/, "_").trim()
                log.warn "Missing target_ref_name for sample ${row.sample_name}; using singlegenome_species (${target_ref_name}) to select the closest target reference."
            }
            target_ref_name = target_ref_name.replaceAll(/\s+/, "_").replaceAll(/\.(fna|fa|fasta)$/, "")

            def accessions = row.archive_accession.split(",")
            return [species_key, row.sample_name, target_ref_name, accessions]
        }

    ch_per_accession = ch_samples
        .flatMap { species, sample, target_ref_name, accessions ->
            accessions.collect { acc -> [species, sample, target_ref_name, acc.trim()] }
        }

    // ── Download ────────────────────────────────────────────────────────────
    DOWNLOAD_FASTQ(ch_per_accession)

    ch_per_sample = DOWNLOAD_FASTQ.out.fastq
        .groupTuple(by: [0, 1, 2])
        .map { species, sample, target_ref_name, fastq_lists ->
            [species, sample, target_ref_name, fastq_lists.flatten()]
        }

    MERGE_FASTQ(ch_per_sample)

    // ── Trim ────────────────────────────────────────────────────────────────
    FASTP(MERGE_FASTQ.out.merged)

    // ── Assembly channels ───────────────────────────────────────────────────
    ch_yer = Channel.fromPath("${params.assemblies_Yersinia}/*.fna")
        .map { f -> ["Yersinia", f.baseName, f] }

    ch_hep = Channel.fromPath("${params.assemblies_Hepatitis}/*.fna")
        .map { f -> ["Hepatitis", f.baseName, f] }

    ch_her = Channel.fromPath("${params.assemblies_Herpesvirus}/*.fna")
        .map { f -> ["Herpesvirus", f.baseName, f] }

    ch_sal = Channel.fromPath("${params.assemblies_Salmonella}/*.fna")
        .map { f -> ["Salmonella", f.baseName, f] }

    ch_par = Channel.fromPath("${params.assemblies_Parvovirus}/*.fna")
        .map { f -> ["Parvovirus", f.baseName, f] }

    ch_myc = Channel.fromPath("${params.assemblies_Mycobacterium}/*.fna")
        .map { f -> ["Mycobacterium", f.baseName, f] }

    ch_var = Channel.fromPath("${params.assemblies_Variola}/*.fna")
        .map { f -> ["Variola", f.baseName, f] }

    ch_assemblies = ch_yer.mix(ch_hep, ch_her, ch_sal, ch_par, ch_myc, ch_var)

    // ── Build indices ───────────────────────────────────────────────────────
    BOWTIE2_INDEX(ch_assemblies)
    BWA_INDEX(ch_assemblies)

    // ── Target-confirmed read extraction, once per aligner ──────────────────
    // The target reference is selected either by exact basename match or,
    // when target_ref_name is inferred from singlegenome_species, by species prefix.
    // Example: target_ref_name = Yersinia_pestis matches
    //          Yersinia_pestis__GCF_000009065.1.fna
    def isTargetRef = { target_ref_name, ref_name ->
        def targetStem = target_ref_name.replaceAll(/\s+/, "_").replaceAll(/\.(fna|fa|fasta)$/, "")
        ref_name == targetStem || ref_name.startsWith(targetStem + "__")
    }

    ch_target_bwa_input = FASTP.out.trimmed
        .combine(BWA_INDEX.out.index, by: 0)
        .filter { species, sample, target_ref_name, reads, ref_name, ref_fasta, bwa_idx ->
            isTargetRef(target_ref_name, ref_name)
        }

    ch_target_bt2_input = FASTP.out.trimmed
        .combine(BOWTIE2_INDEX.out.index, by: 0)
        .filter { species, sample, target_ref_name, reads, ref_name, ref_fasta, bt2_idx ->
            isTargetRef(target_ref_name, ref_name)
        }

    TARGET_BWA_MAP(ch_target_bwa_input)
    TARGET_BOWTIE2_MAP(ch_target_bt2_input)

    // Both target BAMs go through the same filter process; the aligner name
    // travels in the tuple so the two read sets never get mixed up.
    FILTER_FASTQ_TO_TARGET_READS(
        TARGET_BWA_MAP.out.target_bam.mix(TARGET_BOWTIE2_MAP.out.target_bam)
    )

    ch_target_reads = FILTER_FASTQ_TO_TARGET_READS.out.target_reads
        .branch {
            bwa:     it[3] == "BWA"
            bowtie2: it[3] == "Bowtie2"
        }

    // ── Remap each aligner's own target-confirmed reads to all references ───
    ch_bt2_map_input = ch_target_reads.bowtie2
        .combine(BOWTIE2_INDEX.out.index, by: 0)
        .map { species, sample, target_ref_name, aligner, reads, ref_name, ref_fasta, bt2_idx ->
            [species, sample, reads, ref_name, ref_fasta, bt2_idx]
        }

    ch_bwa_map_input = ch_target_reads.bwa
        .combine(BWA_INDEX.out.index, by: 0)
        .map { species, sample, target_ref_name, aligner, reads, ref_name, ref_fasta, bwa_idx ->
            [species, sample, reads, ref_name, ref_fasta, bwa_idx]
        }

    // ── Map ─────────────────────────────────────────────────────────────────
    BOWTIE2_MAP(ch_bt2_map_input)
    BWA_MAP(ch_bwa_map_input)

    ch_all_bams = BOWTIE2_MAP.out.bam.mix(BWA_MAP.out.bam)

    // ── Stats ───────────────────────────────────────────────────────────────
    MAPPING_STATS(ch_all_bams)

    // ── Subsample ───────────────────────────────────────────────────────────
    ch_subsample_input = MAPPING_STATS.out.stats
        .combine(Channel.fromList(params.subsample_depths))

    SUBSAMPLE_BAM(ch_subsample_input)
    PASSTHROUGH_FULL_BAM(MAPPING_STATS.out.stats)

    ch_all_for_dp = SUBSAMPLE_BAM.out.sub_bam
        .mix(PASSTHROUGH_FULL_BAM.out.full_bam)

    // ── DamageProfiler ──────────────────────────────────────────────────────
    DAMAGE_PROFILER(ch_all_for_dp)

    // ── Compute metrics ─────────────────────────────────────────────────────
    ch_metrics_input = DAMAGE_PROFILER.out.profiles
        .map { species, sample, ref, aligner, n, dp5, dp3, stats ->
            [species, sample, ref, aligner, n, dp5, dp3, stats]
        }

    ch_scripts = Channel.fromPath(params.scripts_dir)
    COMPUTE_METRICS(ch_metrics_input, ch_scripts.first())

    // ── Merge metrics, written to results root ──────────────────────────────
    ch_grouped_metrics = COMPUTE_METRICS.out.metrics
        .groupTuple(by: 0)

    MERGE_METRICS(ch_grouped_metrics)
}
