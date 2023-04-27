version 1.0

workflow CountTheBeans {
    input {
        File  bam
        File? bai
        Boolean use_local_ssd
    }

    call Count {input: bam = bam, bai = bai, disk_type = if(use_local_ssd) then "LOCAL" else "SSD"}

    output {
        Map[String, Int] bean_counts = {'raw': Count.raw_count,
                                        'non_2304': Count.non_2304_count,
                                        'mmml': Count.bean_count,
                                        'non_2304_mml': Count.non_2304_bean_count,}
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
