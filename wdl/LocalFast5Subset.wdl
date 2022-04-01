version 1.0
import "tasks/Utils.wdl" as Utils
import "tasks/Structs.wdl"

workflow LocalFast5 {
    input {
        Array[String]+ loci
        String aligned_bam
        File   aligned_bai
        String gcs_output_dir
        File summary_txt
        String prefix
    }

    scatter (locus in loci) {
        call Utils.SubsetBam {
            input:
                bam = aligned_bam,
                bai = aligned_bai,
                locus = locus
        }
    }

    if (length(loci) > 1)  {
        call Utils.MergeBams {
            input:
                bams = SubsetBam.subset_bam,
                prefix = "merged"
        }
    }

    File subset_bam = select_first([MergeBams.merged_bam, SubsetBam.subset_bam[0]])

    call GetReadnames {
        input:
            bam = subset_bam
    }

    call GetFast5Filenames {
        input:
            readnames = GetReadnames.readnames,
            summary_file = summary_txt
    }

    call GetLocalFast5 {
        input:
            readnames = GetReadnames.readnames,
            filenames = GetFast5Filenames.filenames,
            numfiles = GetFast5Filenames.numlines,
            gcs_output_dir = gcs_output_dir,
            prefix = prefix
    }
}

task GetLocalFast5 {
    input {
        File readnames
        File filenames
        Int numfiles
        String gcs_output_dir
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = ceil(0.6 * numfiles)

    command <<<
        set -euxo pipefail
        num_core=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)

        mkdir fast5
        mkdir output

        while read filename; do gsutil cp $filename fast5/ ; done < ~{filenames}

        fast5_subset -i fast5 -s output -l ~{readnames} -t $num_core

        ## save output
        cd output
        gsutil -m cp *.fast5 ~{gcs_output_dir}/~{prefix}/
    >>>

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          8,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-ont:latest"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task GetFast5Filenames {
    input {
        File readnames
        File summary_file

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 3*ceil((size(readnames, "GB")+size(summary_file, "GB")))

    String fast5_dir = sub(summary_file, basename(summary_file), "fast5_pass/")

    command <<<
        set -euxo pipefail

        sort -k3,3 ~{summary_file} -o ~{summary_file}
        join -1 1 -2 3 -o 2.2 ~{readnames} ~{summary_file} |sort -u > filenames
        sed -i '' 's/^/~{fast5_dir}/' filenames
        wc -l filenames | awk '{print $1}'
    >>>

    output {
        File filenames = "filenames"
        Int numlines = read_int(stdout())
    }
    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "ubuntu:16.04"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task GetReadnames {
    input {
        File bam

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail

        samtools view ~{bam} | cut -f1| sort -u > readnames
    >>>

    output {
        File readnames = "readnames"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.9"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}