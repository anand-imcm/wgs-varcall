version 1.0

import "./tasks/deepvariant.wdl" as dvar
import "./tasks/gauchian.wdl" as gvar
import "./tasks/summary.wdl" as stats
import "./tasks/annotation.wdl" as ann

workflow main {
    String pipeline_version = "1.0.0"
    String container_src = "docker.io/library/dev-wgs:~{pipeline_version}"
    String deepvariant_docker = "ghcr.io/anand-imcm/deepvariant:1.9.0"
    String vep_docker = "ghcr.io/anand-imcm/ensembl-vep:release_113.3"

    input {
        File cram
        String sample_id
        File genome_reference
        Boolean use_gauchian = true
        Boolean use_deepvariant = false
        Boolean use_vep = false
        File vep_cache
    }
    if(use_gauchian){
        call gvar.gauchian {
            input:
                cram = cram,
                genome_reference = genome_reference,
                sample_id = sample_id,
                docker = container_src
        }
    }
    if(use_deepvariant){
        call dvar.deepvariant {
            input:
                cram = cram,
                genome_reference = genome_reference,
                sample_id = sample_id,
                docker = deepvariant_docker
        }
        if (use_vep){
            call ann.vep {
                input: 
                    vcf = deepvariant.filtered_vcf,
                    cache = vep_cache,
                    genome_reference = genome_reference,
                    file_label = sample_id,
                    docker = vep_docker
            }
        }
    }
    call stats.summary {
        input:
            cram = cram,
            reference = genome_reference,
            sample_id = sample_id,
            vcf = deepvariant.filtered_vcf,
            vep_stats = vep.stats,
            vep_out = vep.annotated_vcf,
            docker = container_src
    }
    output {
        File? gauchian_summary = gauchian.gauchian_tsv
        File? gauchian_json = gauchian.gauchian_json
        File? deepvariant_all_vcf = deepvariant.all_variants_vcf
        File? deepvariant_all_vcf_stats = deepvariant.all_vcf_stats
        File? deepvariant_summary = deepvariant.summary
        File? deepvariant_gvcf = deepvariant.all_variants_gvcf
        File? deepvariant_filtered_vcf = deepvariant.filtered_vcf
        File? deepvariant_filtered_vcf_stats = summary.vcf_stats
        File? vep_annotated_vcf = vep.annotated_vcf
        File? vep_annotated_tsv = summary.annotation_tsv
        File? vep_annotation_stats = vep.stats
        File qc_report = summary.multiqc_report
        File qc_plots = summary.multiqc_plots
        File qc_data = summary.multiqc_data
        File alignment_summary = summary.cram_summary
        File depth_summary = summary.depth_summary
        File global_depth_summary = summary.global_depth
    }
    meta {
        description: "A WDL workflow for variant calling and annotation from WGS CRAM files. It uses Google's DeepVariant for genome-wide SNV and indel detection, and Illumina's Gauchian for precise variant calling in the GBA gene region"
        author: "Anand Maurya"
        email: "anand.maurya@well.ox.ac.uk"
    }
}