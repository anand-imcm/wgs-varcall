version 1.0

task summary {
    input {
        File cram
        File reference
        String sample_id
        File? vcf
        File? vep_stats
        File? vep_out
        Int memory_gb = 48
        Int cpu = 32
        String docker
    }  

    command <<<
        set -euo pipefail
        ref=$(basename "~{reference}")
        ln -s ~{reference} $ref
        samtools faidx -@ $(( $(nproc) * 3 / 4 )) $ref
        ln -s ~{cram} ~{sample_id}.cram
        samtools index -@ $(( $(nproc) * 3 / 4 )) ~{sample_id}.cram
        if [ -n "~{vcf}" ]; then
            ln -s ~{vcf} ~{sample_id}.vcf.gz
            bcftools stats \
                ~{sample_id}.vcf.gz \
                --fasta-ref $ref > ~{sample_id}.vcf.stats.txt
        fi
        if [ -n "~{vep_stats}" ]; then
            ln -s ~{vep_stats} ~{sample_id}.vep.txt
        fi
        if [ -n "~{vep_out}" ]; then
            headers="chrom\tpos\tref\talt\tqual\tfilter\tgenotype\tgenotype_qual\tread_depth\tallele_depth\tvariant_allele_frac\tgenotype_likelihood\tVEP_Allele\tConsequence\tIMPACT\tSYMBOL\tGene\tFeature_type\tFeature\tBIOTYPE\tEXON\tINTRON\tHGVSc\tHGVSp\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tDISTANCE\tSTRAND\tFLAGS\tSYMBOL_SOURCE\tHGNC_ID\tCANONICAL\tMANE\tMANE_SELECT\tMANE_PLUS_CLINICAL\tCCDS\tENSP\tRefSeq\tREFSEQ_MATCH\tSOURCE\tREFSEQ_OFFSET\tGIVEN_REF\tUSED_REF\tBAM_EDIT\tSIFT\tPolyPhen\tHGVS_OFFSET\tAF\tAFR_AF\tAMR_AF\tEAS_AF\tEUR_AF\tSAS_AF\tgnomADg_AF\tgnomADg_AFR_AF\tgnomADg_AMI_AF\tgnomADg_AMR_AF\tgnomADg_ASJ_AF\tgnomADg_EAS_AF\tgnomADg_FIN_AF\tgnomADg_MID_AF\tgnomADg_NFE_AF\tgnomADg_REMAINING_AF\tgnomADg_SAS_AF\tMAX_AF\tMAX_AF_POPS\tCLIN_SIG\tSOMATIC\tPHENO\tPUBMED\tClinVar\tClinVar_CLNSIG\tClinVar_CLNREVSTAT\tClinVar_CLNDN"
            echo -e $headers > ~{sample_id}.filtered.annotated.variants.tsv
            bcftools +split-vep ~{vep_out} -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t[%GT\t%GQ\t%DP\t%AD\t%VAF\t%PL]\t%CSQ\n' -d -A tab >> ~{sample_id}.filtered.annotated.variants.tsv
        fi
        mosdepth --threads $(( $(nproc) * 3 / 4 )) \
            --no-per-base \
            --fasta $ref \
            ~{sample_id} \
            ~{sample_id}.cram
        samtools flagstat \
            -@ $(( $(nproc) * 3 / 4 )) \
            ~{sample_id}.cram > ~{sample_id}.flagstat.txt
        multiqc . \
            --force \
            --no-ai \
            --filename ~{sample_id} \
            --export \
            --zip-data-dir
        zip -r ~{sample_id}.plots.zip ~{sample_id}_plots/
        mv ~{sample_id}_data.zip ~{sample_id}.data.zip
    >>>

    output {
        File cram_summary = sample_id + ".flagstat.txt"
        File depth_summary = sample_id + ".mosdepth.summary.txt"
        File global_depth = sample_id + ".mosdepth.global.dist.txt"
        File? vcf_stats = sample_id + ".vcf.stats.txt"
        File? annotation_tsv = sample_id + ".filtered.annotated.variants.tsv"
        File multiqc_report = sample_id + ".html"
        File multiqc_data = sample_id + ".data.zip"
        File multiqc_plots = sample_id + ".plots.zip"
    }

    runtime {
        docker: "~{docker}"
        cpu: "~{cpu}"
        memory: "~{memory_gb}GB"
        disks: "local-disk ~{ceil(size([cram, reference], 'GB')) * 3} HDD"
    }
}