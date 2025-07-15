version 1.0

task gauchian {
    input {
        File cram
        File genome_reference
        String sample_id
        String docker
        Int memory_gb = 48
        Int cpu = 32
    }  

    command <<<
        set -euo pipefail
        ref=$(basename "~{genome_reference}")
        ln -s ~{genome_reference} $ref
        samtools faidx -@ $(( $(nproc) * 3 / 4 )) $ref
        ln -s ~{cram} ~{sample_id}.cram
        samtools index -@ $(( $(nproc) * 3 / 4 )) ~{sample_id}.cram
        echo ~{sample_id}.cram > manifest.txt
        gauchian \
            --manifest manifest.txt \
            --genome 38 \
            --reference $ref \
            --prefix ~{sample_id}.gauchian \
            --outDir . \
            --threads $(nproc)
    >>>

    output {
        File gauchian_tsv = sample_id + ".gauchian.tsv"
        File gauchian_json = sample_id + ".gauchian.json"
    }
    runtime {
        docker: "~{docker}"
        bootDiskSizeGb: 20
        cpu: "~{cpu}"
        memory: "~{memory_gb}GB"
        disks: "local-disk ~{ceil(size([cram,genome_reference], 'GB')) * 4} HDD"
    }
}