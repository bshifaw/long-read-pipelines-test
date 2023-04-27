version 1.0

workflow CountTheBeans {
    input {
        File  bam
        File? bai
        Boolean use_local_ssd
    }

    call Count { input: bam = bam, bai = bai, disk_type = if(use_local_ssd) then "LOCAL" else "SSD"}

    call GatherBitter { input: bam = bam, bai = bai, disk_type = if(use_local_ssd) then "LOCAL" else "SSD"}

    output {
        Map[String, Int] bean_counts = {'raw': Count.raw_count,
                                        'non_2304': Count.non_2304_count,
                                        'mmml': Count.bean_count,
                                        'non_2304_mml': Count.non_2304_bean_count,}
        Map[String, File] bitter_reads = {'no_ml': GatherBitter.no_ml_reads,
                                          'no_ml': GatherBitter.no_mm_reads,
                                          'missing_one_tag': GatherBitter.reads_missing_only_one_tag,
                                          'missing_both_tags': GatherBitter.reads_missing_both_tags}
    }
}

task Count {
    input {
        File  bam
        File? bai
        String disk_type
    }

    output {
        Int raw_count  = read_int("raw_count.txt")
        Int non_2304_count = read_int("non_2304_count.txt")
        Int bean_count = read_int("bean_count.txt")
        Int non_2304_bean_count = read_int("non_2304_bean_count.txt")
    }

    command <<<
        set -eux

        samtools view -@1 -c ~{bam} > raw_count.txt &
        samtools view -@1 -c -F 2304 ~{bam} > non_2304_count.txt &

        samtools view -@1 ~{bam} | grep "ML:B:C" | grep -c "MM:Z" > bean_count.txt &
        samtools view -@1 -F 2304  ~{bam} | grep "ML:B:C" | grep -c "MM:Z" > non_2304_bean_count.txt &

        wait
    >>>

    runtime {
        cpu:            10
        memory:         "48 GiB"
        disks:          "local-disk 375 ~{disk_type}"
        preemptible:    2
        maxRetries:     1
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}

task GatherBitter {
    input {
        File  bam
        File? bai
        String disk_type
    }

    String p = basename(bam, ".bam")

    output {
        File no_ml_reads = "~{p}.no_ML.bam"
        File no_mm_reads = "~{p}.no_MM.bam"

        File reads_missing_only_one_tag = "missing_only_one_tag.read_names.txt"
        File reads_missing_both_tags    = "no_mm_and_ml.read_names.txt"
    }

    command <<<
        set -eux

        samtools view -@3 -h ~{bam} \
            | grep -v "ML:B:C" \
        > "~{p}.no_ML.bam" &

        samtools view -@3 -h ~{bam} \
            | grep -v "MM:Z" \
        > "~{p}.no_MM.bam" &

        wait

        samtools view "~{p}.no_ML.bam" | awk -F '\t' '{print $1}' | sort > no_ml.txt
        samtools view "~{p}.no_MM.bam" | awk -F '\t' '{print $1}' | sort > no_mm.txt
        comm -3 \
            no_ml.txt \
            no_mm.txt
        > "missing_only_one_tag.read_names.txt"
        comm -12 \
            no_ml.txt \
            no_mm.txt
        > "no_mm_and_ml.read_names.txt"
    >>>

    runtime {
        cpu:            10
        memory:         "48 GiB"
        disks:          "local-disk 375 ~{disk_type}"
        preemptible:    2
        maxRetries:     1
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }
}
