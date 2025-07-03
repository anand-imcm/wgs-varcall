version 1.0

task summary {
    input {
        File cram
        File reference
        String sample_id
        File? vcf
        File? vep_stats
        Int memory_gb = 24
        Int cpu = 16
        String docker
    }  

    command <<<
        set -euo pipefail
        samtools index ~{cram}
        samtools faidx ~{reference}
        if [ -n "~{vcf}" ]; then
            ln -s ~{vcf} ~{sample_id}.vcf.gz
            bcftools stats \
                ~{sample_id}.vcf.gz \
                --fasta-ref ~{reference} > ~{sample_id}.vcf.stats.txt
        fi
        if [ -n "~{vep_stats}" ]; then
            ln -s ~{vep_stats} ~{sample_id}.vep.txt
        fi
        samtools flagstat \
            ~{cram} > ~{sample_id}.flagstat.txt
        mosdepth \
            --threads $(nproc) \
            --no-per-base \
            --fasta ~{reference} \
            ~{sample_id} \
            ~{cram}
        multiqc . \
            --force \
            --no-ai \
            --filename ~{sample_id} \
            --export \
            --zip-data-dir
        zip -r ~{sample_id}_plots.zip ~{sample_id}_plots/
    >>>

    output {
        File cram_index = cram + ".crai"
        File genome_index = reference + ".fai"
        File cram_summary = sample_id + ".flagstat.txt"
        File depth_summary = sample_id + ".mosdepth.summary.txt"
        File global_depth = sample_id + ".mosdepth.global.dist.txt"
        File? vcf_stats = sample_id + ".vcf.stats.txt"
        File multiqc_report = sample_id + ".html"
        File multiqc_data = sample_id + "_data.zip"
        File multiqc_plots = sample_id + "_plots.zip"
    }

    runtime {
        docker: "~{docker}"
        cpu: "~{cpu}"
        memory: "~{memory_gb}GB"
        disks: "local-disk ~{ceil(size(cram, 'GB')) * 2} HDD"
    }
}