version 1.0

##########################################################################################
## A workflow that performs CCS correction on PacBio HiFi reads from a single flow cell.
## The workflow shards the subreads into clusters and performs CCS in parallel on each cluster.
## Ultimately, all the corrected reads (and uncorrected) are gathered into a single BAM.
## Various metrics are produced along the way.
##########################################################################################

import "tasks/PBUtils.wdl" as PB
import "tasks/Utils.wdl" as Utils
import "tasks/Finalize.wdl" as FF

workflow PBFlowcell {
    input {
        File bam
        File pbi
        File ref_map_file

        String participant_name
        String sample_name

        Int num_shards = 300

        Boolean is_corrected
        Boolean is_ccs
        Boolean is_isoseq

        String gcs_out_root_dir
    }

    parameter_meta {
        bam:                "GCS path to raw subread bam"
        pbi:                "GCS path to pbi index for raw subread bam"
        ref_map_file:       "table indicating reference sequence and auxillary file locations"

        SM:                 "the value to place in the BAM read group's SM field"
        ID:                 "the value to place in the BAM read group's ID field"

        is_corrected:       "indicates whether input has already been CCS-corrected"
        is_ccs:             "indicates whether input represents HiFi data (that might not have been corrected yet)"
        is_isoseq:          "indicates whether input represents IsoSeq data"

        gcs_out_root_dir:   "GCS bucket to store the reads, variants, and metrics files"
    }

    Map[String, String] ref_map = read_map(ref_map_file)
    Map[String, String] map_presets = {
        'CLR': 'SUBREAD',
        'CCS': 'CCS',
        'IsoSeq': 'ISOSEQ'
    }

    # is_corrected  is_ccs  is_isoseq    correct      align
    # 0             0       0            no           SUBREAD
    # 0             1       0            yes          CCS
    # 1             0       0            no           CCS
    # 1             1       0            no           CCS
    # 0             0       1            no           ISOSEQ
    # 0             1       1            yes          ISOSEQ
    # 1             0       1            no           ISOSEQ
    # 1             1       1            no           ISOSEQ

    Boolean correct = !is_corrected && is_ccs
    String experiment_type = "SUBREAD"


    call PB.GetRunInfo { input: bam = bam, SM = sample_name }
    String ID = GetRunInfo.run_info["PU"]

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/PBFlowcell/~{ID}"

    # break one raw BAM into fixed number of shards
    call PB.ShardLongReads { input: unaligned_bam = bam, unaligned_pbi = pbi, num_shards = num_shards }

    # then perform correction on each of the shard
    scatter (subreads in ShardLongReads.unmapped_shards) {
        if (experiment_type != "CLR") {
            call PB.CCS { input: subreads = subreads }
        }

        File unaligned_bam = select_first([CCS.consensus, subreads])

        call PB.Align as AlignReads {
            input:
                bam         = unaligned_bam,
                ref_fasta   = ref_map['fasta'],
                sample_name = participant_name,
                map_preset  = map_presets[experiment_type]
        }
    }

    # merge corrected, unaligned reads
    if (experiment_type != "CLR") {
        call Utils.MergeBams as MergeCCSUnalignedReads { input: bams = select_all(CCS.consensus) }
        call PB.PBIndex as IndexCCSUnalignedReads { input: bam = MergeCCSUnalignedReads.merged_bam }

        call PB.MergeCCSReports as MergeCCSReports { input: reports = select_all(CCS.report), prefix = ID }
        call PB.SummarizeCCSReport { input: report = MergeCCSReports.report }

        call FF.FinalizeToFile as FinalizeCCSUnalignedBam {
            input:
                file    = MergeAlignedReads.merged_bam,
                outfile = outdir + "/reads/ccs/unaligned" + basename(MergeCCSUnalignedReads.merged_bam)
        }

        call FF.FinalizeToFile as FinalizeCCSUnalignedPbi {
            input:
                file    = IndexCCSUnalignedReads.pbi,
                outfile = outdir + "/reads/ccs/unaligned/" + basename(IndexCCSUnalignedReads.pbi)
        }

        call FF.FinalizeToFile as FinalizeCCSReport {
            input:
                file    = MergeCCSReports.report,
                outfile = outdir + "/reads/ccs/unaligned/" + basename(MergeCCSReports.report)
        }
    }

    # merge the corrected per-shard BAM/report into one, corresponding to one raw input BAM
    call Utils.MergeBams as MergeAlignedReads { input: bams = AlignReads.aligned_bam, prefix = ID }
    call PB.PBIndex as IndexAlignedReads { input: bam = MergeAlignedReads.merged_bam }

    call PB.SummarizePBI as SummarizeSubreadsPBI { input: pbi = pbi }
    call PB.SummarizePBI as SummarizeAlignedPBI { input: pbi = IndexAlignedReads.pbi }
    call PB.SummarizePBI as SummarizeAlignedQ20PBI { input: pbi = IndexAlignedReads.pbi, qual_threshold = 20 }

    call Utils.ComputeGenomeLength { input: fasta = ref_map['fasta'] }

    # Finalize data
    String dir_prefix = if (experiment_type != "CLR") then "reads/ccs/aligned" else "reads/subreads/aligned"

    call FF.FinalizeToFile as FinalizeAlignedBam {
        input:
            file    = MergeAlignedReads.merged_bam,
            outfile = outdir + "/" + dir_prefix + "/" + basename(MergeAlignedReads.merged_bam)
    }

    call FF.FinalizeToFile as FinalizeAlignedBai {
        input:
            file    = MergeAlignedReads.merged_bai,
            outfile = outdir + "/" + dir_prefix + "/" + basename(MergeAlignedReads.merged_bai)
    }

    call FF.FinalizeToFile as FinalizeAlignedPbi {
        input:
            file    = IndexAlignedReads.pbi,
            outfile = outdir + "/" + dir_prefix + "/" + basename(IndexAlignedReads.pbi)
    }

    output {
        File? ccs_bam = MergeCCSUnalignedReads.merged_bam
        File? ccs_pbi = IndexCCSUnalignedReads.pbi

        File aligned_bam = MergeAlignedReads.merged_bam
        File aligned_pbi = IndexAlignedReads.pbi

        Float num_records = SummarizeSubreadsPBI.results['reads']
        Float total_bases = SummarizeSubreadsPBI.results['bases']
        Float raw_yield = SummarizeSubreadsPBI.results['yield']
        Float raw_est_fold_cov = SummarizeSubreadsPBI.results['yield']/ComputeGenomeLength.length

        Float polymerase_mean = SummarizeSubreadsPBI.results['polymerase_mean']
        Float polymerase_n50 = SummarizeSubreadsPBI.results['polymerase_n50']

        Float subread_mean = SummarizeSubreadsPBI.results['subread_mean']
        Float subread_n50 = SummarizeSubreadsPBI.results['subread_n50']

        Float ccs_num_records = SummarizeAlignedPBI.results['reads']
        Float ccs_total_length = SummarizeAlignedPBI.results['bases']
        Float ccs_mean_qual = SummarizeAlignedPBI.results['mean_qual']
        Float ccs_yield = SummarizeAlignedPBI.results['yield']
        Float ccs_est_fold_cov = SummarizeAlignedPBI.results['yield']/ComputeGenomeLength.length

        Float ccs_num_records_q20 = SummarizeAlignedQ20PBI.results['reads']
        Float ccs_total_length_q20 = SummarizeAlignedQ20PBI.results['bases']
        Float ccs_mean_qual_q20 = SummarizeAlignedQ20PBI.results['mean_qual']
        Float ccs_yield_q20 = SummarizeAlignedQ20PBI.results['yield']

        File? ccs_report = MergeCCSReports.report
        Float? zmws_input = SummarizeCCSReport.zmws_input
        Float? zmws_pass_filters = SummarizeCCSReport.zmws_pass_filters
        Float? zmws_fail_filters = SummarizeCCSReport.zmws_fail_filters
        Float? zmws_pass_filters_pct = SummarizeCCSReport.zmws_pass_filters_pct
        Float? zmws_fail_filters_pct = SummarizeCCSReport.zmws_fail_filters_pct
    }
}
