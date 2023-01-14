version 1.0

import "tasks/Utils.wdl"
import "tasks/utils/BAMutils.wdl" as BU
import "tasks/utils/GeneralUtils.wdl" as GU
import "tasks/Finalize.wdl" as FF

workflow SplitMergedBamByReadGroup {
    meta {
        description: "Split a BAM file that was aggregated, for the same sample, into pieces by read group."
    }

    input {
        File input_bam
        File? input_bai

        String platform
        Boolean unmap_bam
        Boolean convert_to_fq = false

        String gcs_out_root_dir

        Boolean debug_mode = false
    }

    parameter_meta {
        input_bam: "BAM to be split by read group; doesn't necessarily need to be aligned or sorted."
        input_bai: "(optional) BAI accompanying the BAM"
        platform:  "long reads platform the BAM was generated on; must be one of [PB, ONT]"
        convert_to_fq: "user option to convert to FASTQ (gz) or not"
        gcs_out_root_dir: "place to store the result files"
    }

    if (platform!="PB" && platform!="ONT") {
        String formatted_reason = "Provided value for 'platform' (" + platform + ") isn't supported. Must be one of [PB, ONT]."
        call Utils.StopWorkflow { input: reason = formatted_reason }
    }

    String workflow_name = "SplitMergedBamByReadGroup"
    String outdir = sub(gcs_out_root_dir, "/$", "") + "/~{workflow_name}/" + basename(input_bam, '.bam')

    call BU.GatherBamMetadata {
        input: bam = input_bam
    }

    if (debug_mode) {
        call HowManyReadGroupsAreInUse {input: bam = input_bam}
    }

    # split
    String output_prefix = basename(input_bam, ".bam")
    call Utils.ComputeAllowedLocalSSD as Guess {
        input: intended_gb = 10 + 3*ceil(size([input_bam], "GB"))
    }
    call BU.SplitByRG {
        input:
            bam = input_bam, out_prefix = output_prefix, num_ssds = Guess.numb_of_local_ssd
    }

    scatter (bam in SplitByRG.split_bam) {

        call BU.GetReadGroupInfo { input: uBAM = bam, keys = ['ID', 'LB'], null_value_representation = 'None' }
        String rgid = GetReadGroupInfo.read_group_info['ID']
        String library = GetReadGroupInfo.read_group_info['LB']

        if (debug_mode) {
            call Utils.CountBamRecords { input: bam = bam }
        }

        # drop alignment if so requested
        if (unmap_bam) {
            call BU.UnAlignBam { input: bam = bam, out_prefix = basename(bam, ".bam") }
            call FF.FinalizeToFile as SaveUBam {
                input: file = UnAlignBam.uBAM, outdir = outdir
            }
            if (convert_to_fq) { # save FASTQ only if requested
                call FF.FinalizeToFile as SaveFq {
                    input: file = UnAlignBam.fq, outdir = outdir
                }
            }
            Boolean uBAM_is_empty = UnAlignBam.is_bam_empty
        }
        if (!unmap_bam) {
            call FF.FinalizeToFile as SaveAlnBam {
                input: file = bam, outdir = outdir
            }
        }

        # convert to FASTQ if so requested
        if (convert_to_fq) {
            if (!unmap_bam) { # if alignment unrolling is requested at the same time, it's done there
                call Utils.BamToFastq { input: bam = bam, prefix = basename(bam, ".bam") }
                call FF.FinalizeToFile as SaveFqE {
                    input: file = BamToFastq.reads_fq, outdir = outdir
                }
                if (debug_mode) { call Utils.CountFastqRecords { input: fastq = BamToFastq.reads_fq } }
            }
        }
    }
    Array[String]  phased_rg_ids   = rgid
    Array[String]  phased_bams     = select_first([select_all(SaveUBam.gcs_path), select_all(SaveAlnBam.gcs_path)])
    Array[Boolean?] are_ubams_empty = uBAM_is_empty
    Array[String]? phased_fastqs   = select_first([select_all(SaveFq.gcs_path), select_all(SaveFqE.gcs_path)])

    call GU.CoerceArrayOfPairsToMap as MapRgid2Bams { input: keys = phased_rg_ids, values = phased_bams }
    if (convert_to_fq) {
        call GU.CoerceArrayOfPairsToMap as MapRgid2Fqs { input: keys = phased_rg_ids, values = select_first([phased_fastqs]) }
    }
    if (unmap_bam) {
        call GU.CoerceArrayOfPairsToMap as MapRgid2BamEmptiness { input: keys = phased_rg_ids, values = select_all(are_ubams_empty) }
    }

    call GU.GetTodayDate as today {}

    output {
        Map[String, String] rgid_2_bam = MapRgid2Bams.output_map
        Map[String, Boolean]? rgid_2_ubam_emptyness = MapRgid2BamEmptiness.output_map
        Boolean rgid_2_bam_are_aligned = ! unmap_bam
        Map[String, String]? rgid_2_fastq = MapRgid2Fqs.output_map

        String last_postprocessing_date = today.yyyy_mm_dd
    }
}

task HowManyReadGroupsAreInUse {
    meta {
        desciption: "For debug purposes. Parse each record and collect counts of each read group."
    }
    input {
        File bam
    }

    command <<<
        set -eux

        time \
        samtools view -@3 ~{bam} | \
            grep -Eo "RG:Z:[A-Za-z0-9_/-]+$(printf '\t')" \
        > all_read_readgroup_ids.txt

        time \
        sort all_read_readgroup_ids.txt | \
            uniq -c | \
            awk '{print $2,$1}' | awk '{$1=$1};1' | tr -s ' ' '\t' | sed 's/RG:Z://' \
        > readgroup_ids_and_counts.tsv
    >>>

    output {
        Map[String, Int] rg_ids_freq = read_map("readgroup_ids_and_counts.tsv")
    }

    runtime {
        cpu:    4
        disks:  "local-disk 375 LOCAL"
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}
