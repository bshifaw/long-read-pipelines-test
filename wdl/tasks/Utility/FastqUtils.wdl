version 1.0

import "../../structs/Structs.wdl"


task Stats {
    meta {
        desription:
        "seqkit stats command"
    }
    parameter_meta {
        fastq: "file to collect stats on"
        seq_type: "argument to the --seq-type paramter"
    }

    input {
        File fastq
        String seq_type
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 100 + ceil(size(fastq, "GB"))

    command <<<
        set -eux

        seqkit stats \
            -aT \
            -t ~{seq_type} \
            -o 2_line.tsv ~{fastq}

        datamash transpose \
            < 2_line.tsv \
            | grep -vw "^file" \
            | grep -vw "^format" \
            | grep -vw "^type" \
        > 2.col.map.tsv
    >>>

    output {
        Map[String, Float] res = read_map("2.col.map.tsv")
    }

    ###################
    RuntimeAttr default_attr = object {
        cpu_cores:             4,
        mem_gb:                16,
        disk_gb:               disk_size,
        boot_disk_gb:          10,
        preemptible_tries:     0,
        max_retries:           0,
        docker:                "us.gcr.io/broad-dsp-lrma/lr-seqkit:2.4.0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                   select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory:                select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:        select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible:           select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:            select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker:                select_first([runtime_attr.docker, default_attr.docker])
    }
}

task FilterByLen {
    meta {
        desciption:
        ""
    }
    parameter_meta {
        exclude_len_threshold: "Sequeces shorter than this will be dropped from analysis."
    }

    input {
        File fastq
        Int exclude_len_threshold
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10 + 2*ceil(size(fastq, "GB"))

    String base = basename(basename(fastq, ".fastq.gz"), ".fq.gz")
    String out_prefx = base + ".RL_gt_" + exclude_len_threshold

    command <<<
        set -eux

        seqtk seq \
            -L ~{exclude_len_threshold} \
            ~{fastq} \
        | gzip \
        > "~{out_prefx}.fastq.gz"
    >>>

    output {
        File res = "~{out_prefx}.fastq.gz"
    }

    ###################
    RuntimeAttr default_attr = object {
        cpu_cores:             4,
        mem_gb:                16,
        disk_gb:               disk_size,
        boot_disk_gb:          10,
        preemptible_tries:     0,
        max_retries:           0,
        docker:                "us.gcr.io/broad-dsp-lrma/lr-seqtk:1.3"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                   select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory:                select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:        select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible:           select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:            select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker:                select_first([runtime_attr.docker, default_attr.docker])
    }
}
