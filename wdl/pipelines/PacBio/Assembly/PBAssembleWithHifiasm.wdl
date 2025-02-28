version 1.0

import "../../../tasks/Utility/Utils.wdl" as Utils
import "../../../tasks/Assembly/Hifiasm.wdl" as HA
import "../../../tasks/QC/Quast.wdl" as QuastEval
import "../../../tasks/Utility/Finalize.wdl" as FF

workflow PBAssembleWithHifiasm {

    meta {
        description: "A workflow that performs single sample genome assembly on PacBio HiFi reads from one or more SMRT cells. The multiple SMRT cells data are merged prior to assembly."
    }
    parameter_meta {
        ccs_fqs:            "GCS path to CCS fastq files"

        participant_name:   "name of the participant from whom these samples were obtained"
        prefix:             "prefix for output files"

        ref_fasta_for_eval: "Reference Fasta used for evaluating "
        gcs_out_root_dir:   "GCS bucket to store the reads, variants, and metrics files"
    }

    input {
        Array[File] ccs_fqs

        String participant_name
        String prefix

        File? ref_fasta_for_eval

        String gcs_out_root_dir
    }

    #########################################################################################
    if (length(ccs_fqs) > 1) {
        call Utils.MergeFastqs as MergeAllFastqs { input: fastqs = ccs_fqs }
    }
    File ccs_fq  = select_first([ MergeAllFastqs.merged_fastq, ccs_fqs[0] ])

    call HA.Hifiasm {
        input:
            reads = ccs_fq,
            prefix = prefix
    }

    # todo: assumes ploidy 2
    call QuastEval.Quast as primary_h0_h1_quast {
        input:
            ref = ref_fasta_for_eval,
            is_large = true,
            assemblies = [Hifiasm.primary_tigs,
                          Hifiasm.phased_tigs[0],
                          Hifiasm.phased_tigs[1]]
    }

    call QuastEval.SummarizeQuastReport as primary_h0_h1_quast_summary {
        input: quast_report_txt = primary_h0_h1_quast.report_txt
    }

    #########################################################################################
    # Finalize data
    String workflow_name = "PBAssembleWithHifiasm"

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/" + workflow_name + "/~{prefix}"
    String dir = outdir + "/assembly"

    # assembly results themselves
    call FF.CompressAndFinalize as FinalizeHifiasmPrimaryGFA   { input: outdir = dir, file = Hifiasm.primary_gfa }
    call FF.CompressAndFinalize as FinalizeHifiasmPrimaryFA    { input: outdir = dir, file = Hifiasm.primary_tigs }

    call FF.CompressAndFinalize as FinalizeHifiasmAlternateGFA   { input: outdir = dir, file = Hifiasm.alternate_gfa }
    call FF.CompressAndFinalize as FinalizeHifiasmAlternateFA    { input: outdir = dir, file = Hifiasm.alternate_tigs }

    call FF.FinalizeAndCompress as FinalizeHifiasmHapGFAs  { input: outdir = dir, files = Hifiasm.phased_gfas, prefix = prefix + ".haploGFAs" }
    call FF.FinalizeAndCompress as FinalizeHifiasmHapFAs   { input: outdir = dir, files = Hifiasm.phased_tigs, prefix = prefix + ".haploTigs" }

    call FF.FinalizeToFile as FinalizeQuastReportHtml {
        input: outdir = dir, file = primary_h0_h1_quast.report_html
    }
    call FF.FinalizeAndCompress as FinalizeQuastReports {
        input: outdir = dir, files = primary_h0_h1_quast.report_in_various_formats, prefix = prefix + ".quast_reports"
    }
    call FF.FinalizeToFile as FinalizeQuastSummaryAll {
        input: outdir = dir, file = select_first([primary_h0_h1_quast_summary.quast_metrics_together])
    }
    scatter (report in select_first([primary_h0_h1_quast_summary.quast_metrics]) ) {
        call FF.FinalizeToFile as FinalizeQuastIndividualSummary  { input: outdir = dir, file = report }
    }

    output {
        File hifiasm_primary_gfa  = FinalizeHifiasmPrimaryGFA.gcs_path
        File hifiasm_primary_tigs = FinalizeHifiasmPrimaryFA.gcs_path

        File hifiasm_haploGFAs = FinalizeHifiasmHapGFAs.gcs_path
        File hifiasm_haplotigs = FinalizeHifiasmHapFAs.gcs_path

        File hifiasm_alternate_gfa  = FinalizeHifiasmAlternateGFA.gcs_path
        File hifiasm_alternate_tigs = FinalizeHifiasmAlternateFA.gcs_path

        File? quast_report_html = FinalizeQuastReportHtml.gcs_path
        File? quast_report_in_various_formats = FinalizeQuastReports.gcs_path

        File? quast_summary_on_all = FinalizeQuastSummaryAll.gcs_path

        File? quast_summary_on_primary = FinalizeQuastIndividualSummary.gcs_path[0]
        File? quast_summary_on_H0 = FinalizeQuastIndividualSummary.gcs_path[1]
        File? quast_summary_on_H1 = FinalizeQuastIndividualSummary.gcs_path[2]
    }
}
