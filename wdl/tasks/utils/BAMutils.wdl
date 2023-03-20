version 1.0

import "../Structs.wdl"

task GetReadGroupInfo {
    meta {
        desciption:
        "Get some read group information given a single-readgroup BAM. If the requested keys are absent, a null value is assigned in the returned entry. If the BAM contains multiple read groups, results are undetermined."
    }

    input {
        String uBAM  # not using file as call-caching brings not much benefit

        Array[String] keys
        String null_value_representation = "None"
    }

    parameter_meta {
        keys: "A list of requested fields in the RG line, e.g. ID, SM, LB."
    }

    command <<<
        set -eux

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        samtools view -H ~{uBAM} | grep "^@RG" > one_rg_per_line.txt
        num_rgs=$(wc -l one_rg_per_line.txt | awk '{pritn $1}')
        if [[ num_rgs -gt 1 ]]; then exit 1; fi

        cat one_rg_per_line.txt | tr '\t' '\n' > rh_header.txt

        for attribute in ~{sep=' ' keys}; do
            if grep -q "^${attribute}" rh_header.txt; then
                value=$(grep "^${attribute}" rh_header.txt | awk -F ':' '{print $2}')
            else
                value="~{null_value_representation}"
            fi
            echo -e "${attribute}\t${value}" >> "result.txt"
        done
    >>>

    output {
        Map[String, String] read_group_info = read_map("result.txt")
    }

    runtime {
        cpu:            1
        memory:         "4 GiB"
        disks:          "local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible:    2
        maxRetries:     1
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}

task SplitByRG {
    meta {
        description: "Split a BAM file that was aggregated, for the same sample, into pieces by read group."
    }
    input {
        File bam

        String out_prefix

        Int? num_ssds

        Boolean retain_rgless_records = false
        Boolean sort_and_index = false
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        bam: "BAM to be split"
        out_prefix: "prefix for output bam and bai file names"
        sort_and_index: "if the user wants to (pos-)sort and index the resulting BAMs; this indicates the input BAM is mapped"

        split_bam: "the resuling BAMs, each having reads only in a single read group"
        split_bai: "the accompanying BAIs, if possible and explicit requested"
    }

    Int disk_size = if defined(num_ssds) then 375*select_first([num_ssds]) else 1+3*ceil(size([bam], "GB"))

    Array[String] extra_args = if (retain_rgless_records) then ["-u", "~{out_prefix}_noRG.bam"] else [""]
    command <<<
        set -eux

        samtools view -H ~{bam} | grep "^@RG" > "read_groups_header.txt"
        cat "read_groups_header.txt" | tr '\t' '\n' | grep "^ID:"  | awk -F ':' '{print $2}' > "RG_ids.txt"

        samtools split -@3 \
            -f "~{out_prefix}_%#.bam" \
            ~{sep=" " extra_args} \
            ~{bam}
        if ~{sort_and_index} ;
        then
            # cleanup space for the sorting
            rm ~{bam}
            for split_bam in "~{out_prefix}_"*.bam;
            do
                mv "${split_bam}" temp.bam
                samtools sort \
                    --write-index \
                    -o "${split_bam}##idx##${split_bam}.bai" \
                    temp.bam
            done
        fi
    >>>

    output {
        File read_group_header = "read_groups_header.txt"
        Array[String] rg_ids   = read_lines("RG_ids.txt")
        Array[File]  split_bam = glob("*.bam")
        Array[File?] split_bai = glob("*.bai")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " LOCAL"
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task GatherBamMetadata {
    meta {
        description: "Check several metadata of an input BAM (aliged? sort order? etc)"
    }
    input {
        File bam
    }
    parameter_meta {
        bam: {
            desciption: "BAM to be checked",
            localization_optional: true
        }
    }

    command <<<
        set -eux

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        samtools view -H ~{bam} > header.txt

        grep -F "@HD" header.txt | tr '\t' '\n' > hd.lines.txt
        if grep -q "SO" hd.lines.txt;
        then
            echo "true" > "is_sorted.txt"
            grep "SO" hd.lines.txt | awk -F ':' '{print $2}' > "sort_order.txt"
        else
            echo "false" > "is_sorted.txt"
            echo "NA" > "sort_order.txt"
        fi

        # we use two conditions: @SQ lines in header, and at least some mapped reads
        if grep -q "@SQ" "header.txt";
        then
            export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
            if [[ $(samtools view -F 4 ~{bam} | head | wc -l | awk '{print $1}') -ge 1 ]];
            then
                echo "true" > "is_mapped.txt"
            else
                echo "unsure" > "is_mapped.txt" # this will trigger error in WDL later, but that's intentional because we cannot be sure
            fi
        else
            echo "false" > "is_mapped.txt"
        fi
    >>>

    output {
        Boolean is_aligned = read_boolean("is_mapped.txt")

        Boolean is_sorted = read_boolean("is_sorted.txt")
        String sort_order = read_string("sort_order.txt")
    }

    runtime {
        cpu:    2
        memory: "8 GiB"
        disks:  "local-disk 20 HDD"
        docker: "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.8"
    }
}

task UnAlignBam {
    meta {
        desciption:
        "Go from an aligned BAM to unaligned BAM"
    }
    input {
        File bam

        String out_prefix

        Int? num_ssds
        RuntimeAttr? runtime_attr_override
    }
    Int disk_size = if defined(num_ssds) then 375*select_first([num_ssds]) else 1+3*ceil(size([bam], "GB"))

    String bam_emptyness = "bam_is_empty.txt"

    command <<<
        set -eux

        samtools fastq ~{bam} | gzip > ~{out_prefix}.fq.gz

        # need to reheader to keep history accumulated in the original bam, except those reference contigs
        samtools view -H ~{bam} | grep -v "@SQ" > original_header.txt
        cat original_header.txt

        cnt=$(samtools view -c ~{bam})
        if [[ ${cnt} -eq 0 ]]; then # if input bam is empty, then just output header and signal in the output
            touch "~{out_prefix}.unaligned.bam"
            samtools view -h -o "~{out_prefix}.unaligned.bam" original_header.txt
            echo "true" > ~{bam_emptyness}
            exit 0
        fi
        time \
        samtools import \
            -0 ~{out_prefix}.fq.gz \
            -o ~{out_prefix}.unaligned.tbrh.bam

        time \
        samtools reheader \
            original_header.txt \
            ~{out_prefix}.unaligned.tbrh.bam \
            > ~{out_prefix}.unaligned.bam
        echo "false" > ~{bam_emptyness}
    >>>

    output {
        File fq   = "~{out_prefix}.fq.gz"
        File uBAM = "~{out_prefix}.unaligned.bam"
        Boolean is_bam_empty = read_boolean(bam_emptyness)
    }
    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " LOCAL"
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task ValidateSamFile {
    meta {
        desciption: "Call GATK/Picard ValidateSamFile to validate input BAM: https://bit.ly/3JMutxp."
    }
    parameter_meta {
        validation_mode: "Desired valiation mode"
        disk_type: "Type of disk to use for the computation."
    }

    input {
        File bam
        String validation_mode = "SUMMARY"

        Array[String] validation_errs_to_ignore = ["INVALID_TAG_NM",  # for the purpose we currently have, NM and CIGAR don't matter, and longreads have no mates
                                                    "MISSING_TAG_NM",
                                                    "INVALID_CIGAR",
                                                    "ADJACENT_INDEL_IN_CIGAR",
                                                    "CIGAR_MAPS_OFF_REFERENCE",
                                                    "MISMATCH_MATE_CIGAR_STRING",
                                                    "MATE_CIGAR_STRING_INVALID_PRESENCE",
                                                    "MATE_NOT_FOUND",
                                                    "INVALID_MAPPING_QUALITY",
                                                    "INVALID_FLAG_MATE_UNMAPPED",
                                                    "MISMATCH_FLAG_MATE_UNMAPPED",
                                                    "INVALID_FLAG_MATE_NEG_STRAND",
                                                    "MISMATCH_FLAG_MATE_NEG_STRAND",
                                                    "INVALID_MATE_REF_INDEX",
                                                    "MISMATCH_MATE_REF_INDEX",
                                                    "MISMATCH_MATE_ALIGNMENT_START",
                                                    "MATE_FIELD_MISMATCH",
                                                    "PAIRED_READ_NOT_MARKED_AS_FIRST_OR_SECOND"
                                                   ]

        String disk_type = "LOCAL"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = ceil(size(bam, "GiB")) + 50
    String output_basename = basename(basename(bam, ".bam"), ".cram")
    String output_name = "${output_basename}_${validation_mode}.txt"

    command <<<
        set -eux

        gatk ValidateSamFile \
            --INPUT ~{bam} \
            --OUTPUT ~{output_name} \
            --MODE ~{validation_mode} \
            ~{true="--IGNORE " false="" 0<length(validation_errs_to_ignore)} \
            ~{sep=" --IGNORE " validation_errs_to_ignore}

    >>>

    output {
        File validation_report = "${output_name}"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-gatk/gatk:4.4.0.0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " ~{disk_type}"
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task CountReadGroups {
    meta {
        desciption: "Count the number of RG lines in the header of the BAM file."
    }
    input {
        String bam  # not using file as call-caching brings not much benefit
    }

    command <<<
        set -eux

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        samtools view -H ~{bam} | grep -c "^@RG" > "rg_cnt.txt"
    >>>

    output {
        Int num_rg = read_int("rg_cnt.txt")
    }
    runtime {
        cpu:            1
        memory:         "4 GiB"
        disks:          "local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible:    2
        maxRetries:     1
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}
