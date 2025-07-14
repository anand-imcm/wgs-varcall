version 1.0

task deepvariant {
    input {
        File cram
        File genome_reference
        String sample_id
        String docker = "google/deepvariant:1.9.0"
        Int memory_gb = 48
        Int cpu = 32
    }  

    command <<<
        set -euo pipefail
        export TF_ENABLE_ONEDNN_OPTS=0
        export TF_CPP_MIN_LOG_LEVEL=2
        ref=$(basename "~{genome_reference}")
        ln -s ~{genome_reference} $ref
        samtools faidx $ref
        ln -s ~{cram} ~{sample_id}.cram
        samtools index -@ $(( $(nproc) * 3 / 4 )) ~{sample_id}.cram
        /opt/deepvariant/bin/run_deepvariant \
            --model_type WGS \
            --vcf_stats_report=true \
            --call_variants_extra_args="allow_empty_examples=true" \
            --num_shards $(( $(nproc) * 3 / 4 )) \
            --ref $ref \
            --reads ~{sample_id}.cram \
            --output_vcf ~{sample_id}.vcf.gz \
            --output_gvcf ~{sample_id}.g.vcf.gz
        bcftools filter \
            -i 'FILTER="PASS" && FORMAT/GQ > 20 && FORMAT/DP > 5 && QUAL > 30' \
            ~{sample_id}.vcf.gz \
            -Oz -o ~{sample_id}.pass.filtered.vcf.gz
        bcftools norm \
            ~{sample_id}.pass.filtered.vcf.gz \
            -f $ref \
            -m -any \
            -Oz -o ~{sample_id}.pass.filtered.norm.vcf.gz
        bcftools stats \
            ~{sample_id}.vcf.gz \
            --fasta-ref $ref > ~{sample_id}.all.vcf.stats.txt
    >>>

    output {
        File all_variants_vcf = sample_id + ".vcf.gz"
        File all_variants_gvcf = sample_id + ".g.vcf.gz"
        File summary = sample_id + ".visual_report.html"
        File all_vcf_stats = sample_id + ".all.vcf.stats.txt"
        File filtered_vcf = sample_id + ".pass.filtered.norm.vcf.gz"
    }

    runtime {
        docker: "~{docker}"
        bootDiskSizeGb: 20
        cpu: "~{cpu}"
        memory: "~{memory_gb}GB"
        disks: "local-disk ~{ceil(size([cram,genome_reference], 'GB')) * 3} HDD"
    }
}