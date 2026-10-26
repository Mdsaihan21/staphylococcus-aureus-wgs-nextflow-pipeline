nextflow.enable.dsl=2

params.run_accession = 'ERR16986790'

process DOWNLOAD_READS {
    tag { accession }

    input:
    val accession

    output:
    tuple val(accession), path("${accession}_1.fastq.gz"), path("${accession}_2.fastq.gz")

    script:
    """
    set -euo pipefail

    curl -L --fail --retry 3 \
      "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${accession}&result=read_run&fields=fastq_ftp&format=tsv" \
      -o ena_report.tsv

    awk -F '\\t' 'NR==2 {print \$2}' ena_report.tsv | tr ';' '\\n' > fastq_paths.txt

    test "\$(wc -l < fastq_paths.txt)" -eq 2 || {
      echo "ERROR: Expected two paired-end FASTQ paths"
      exit 1
    }

    i=1
    while IFS= read -r ftp_path; do
      curl -L --fail --retry 3 "https://\${ftp_path}" \
        -o "${accession}_\${i}.fastq.gz"
      i=\$((i + 1))
    done < fastq_paths.txt

    gzip -t "${accession}_1.fastq.gz"
    gzip -t "${accession}_2.fastq.gz"
    """
}

process FASTQC_RAW {
    tag { sample }
    cpus 2

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    path "*_fastqc.html", emit: html
    path "*_fastqc.zip", emit: zip

    script:
    """
    fastqc --threads ${task.cpus} "${read1}" "${read2}"
    """
}

process FASTP_TRIM {
    tag { sample }
    cpus 2

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    tuple val(sample), path("${sample}_1.trimmed.fastq.gz"), path("${sample}_2.trimmed.fastq.gz"), emit: reads
    path "${sample}_fastp.html", emit: html
    path "${sample}_fastp.json", emit: json

    script:
    """
    fastp \
      --in1 "${read1}" \
      --in2 "${read2}" \
      --out1 "${sample}_1.trimmed.fastq.gz" \
      --out2 "${sample}_2.trimmed.fastq.gz" \
      --html "${sample}_fastp.html" \
      --json "${sample}_fastp.json" \
      --thread ${task.cpus} \
      --detect_adapter_for_pe
    """
}

process FASTQC_TRIMMED {
    tag { sample }
    cpus 2

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    path "*_fastqc.html", emit: html
    path "*_fastqc.zip", emit: zip

    script:
    """
    fastqc --threads ${task.cpus} "${read1}" "${read2}"
    """
}

process SPADES_ASSEMBLY {
    tag { sample }
    cpus 1
    memory '2 GB'
    publishDir "${projectDir}/results/assembly", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    tuple val(sample), path("spades_out/contigs.fasta"), emit: contigs
    path "spades_out/scaffolds.fasta", optional: true, emit: scaffolds
    path "spades_out/spades.log", optional: true, emit: log

    script:
    """
    spades.py \
      --pe1-1 "${read1}" \
      --pe1-2 "${read2}" \
      --threads ${task.cpus} \
      --memory 2 \
      --isolate \
      --only-assembler \
      -o spades_out
    """
}

process QUAST_ASSEMBLY {
    tag { sample }
    cpus 1

    publishDir "${projectDir}/results/quast", mode: 'copy'

    input:
    tuple val(sample), path(contigs)

    output:
    path "quast_out", emit: report

    script:
    """
    quast.py "${contigs}" \
      -o quast_out \
      -t ${task.cpus} \
      --min-contig 200
    """
}

process BUSCO_COMPLETENESS {
    tag { sample }
    cpus 1
    memory '3 GB'

    publishDir "${projectDir}/results/busco", mode: 'copy'

    input:
    tuple val(sample), path(contigs)

    output:
    path "busco_out", emit: report

    script:
    """
    busco \
      -i "${contigs}" \
      -m genome \
      --auto-lineage-prok \
      -c ${task.cpus} \
      -o "${sample}" \
      --out_path busco_out
    """
}

process PROKKA_ANNOTATION {
    tag { sample }
    cpus 2
    memory '3 GB'

    publishDir "${projectDir}/results/prokka", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(contigs)

    output:
    path "${sample}.gbk", emit: gbk
    path "${sample}.faa", emit: proteins
    path "${sample}.gff", emit: gff
    path "${sample}.tsv", emit: tsv

    script:
    """
    prokka \
      --outdir prokka_out \
      --prefix ${sample} \
      --cpus ${task.cpus} \
      --kingdom Bacteria \
      --genus Staphylococcus \
      --species aureus \
      --strain ${sample} \
      "${contigs}"

    cp prokka_out/${sample}.gbk .
    cp prokka_out/${sample}.faa .
    cp prokka_out/${sample}.gff .
    cp prokka_out/${sample}.tsv .
    """
}


workflow {
    reads_ch = DOWNLOAD_READS(Channel.value(params.run_accession))
    FASTQC_RAW(reads_ch)
    trimmed_ch = FASTP_TRIM(reads_ch)
    FASTQC_TRIMMED(trimmed_ch.reads)
    SPADES_ASSEMBLY(trimmed_ch.reads)
    QUAST_ASSEMBLY(SPADES_ASSEMBLY.out.contigs)
    BUSCO_COMPLETENESS(SPADES_ASSEMBLY.out.contigs)
    PROKKA_ANNOTATION(SPADES_ASSEMBLY.out.contigs)
}
