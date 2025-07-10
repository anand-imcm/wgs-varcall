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
        samtools index ~{cram}
        samtools faidx ~{genome_reference}
        ln -s ~{cram} ~{sample_id}.cram
        ln -s ~{cram}.crai ~{sample_id}.crai
        echo ~{sample_id}.cram > manifest.txt
        gauchian \
            --manifest manifest.txt \
            --genome 38 \
            --reference ~{genome_reference} \
            --prefix ~{sample_id}.gauchian \
            --outDir . \
            --threads $(nproc)
    >>>

    output {
        File gauchian_tsv = sample_id + ".gauchian.json"
        File gauchian_json = sample_id + ".gauchian.tsv"
    }
    runtime {
        docker: "~{docker}"
        bootDiskSizeGb: 20
        cpu: "~{cpu}"
        memory: "~{memory_gb}GB"
        disks: "local-disk ~{ceil(size([cram,genome_reference], 'GB')) * 3} HDD"
    }
}