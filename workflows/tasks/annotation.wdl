version 1.0

task vep {
    
    input {
        File vcf
        File cache
        File genome_reference
        String file_label
        String docker = "ensemblorg/ensembl-vep:release_110.1"
        Int memory_gb = 48
        Int cpu = 32
    }

    command <<<
        set -euo pipefail
        ln -s ~{genome_reference} genome_reference.fasta
        unzip ~{cache} -d vep_cache/
        perl /opt/vep/src/ensembl-vep/vep --force_overwrite \
            --input_file ~{vcf} \
            --vcf \
            --output_file ~{file_label}.annotated.vcf \
            --stats_text \
            --stats_file ~{file_label}.vep.txt \
            --cache \
            --dir_cache vep_cache/ \
            --merged \
            --fasta genome_reference.fasta \
            --fork $(( $(nproc) * 3 / 4 )) \
            --numbers --offline --hgvs --shift_hgvs 0 --terms SO --symbol \
            --sift b --polyphen b --total_length --ccds --canonical --biotype \
            --protein --xref_refseq --mane --pubmed --af --max_af --af_1kg --af_gnomadg \
            --custom file=vep_cache/clinvar.vcf.gz,short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN
    >>>

    output {
        File annotated_vcf = file_label + ".annotated.vcf"
        File stats = file_label + ".vep.txt"
    }

    runtime {
        docker: "~{docker}"
        cpu: "~{cpu}"
        memory: "~{memory_gb}GB"
        disks: "local-disk ~{ceil(size([cache, genome_reference, vcf], 'GB')) * 3} HDD"
    }
}