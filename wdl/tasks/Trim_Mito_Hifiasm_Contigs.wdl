version 1.0

import "tasks/Structs.wdl"

workflow Trim_Contigs {
    input {
        File assembly_fasta
    }

    call Filter_Contigs {
        input:
        assembly_fasta = assembly_fasta
    }

    call Self_Align {
        input:
        filtered_contigs = Filter_Contigs.filtered_contigs
    }

    call Align_to_Candidates {
        input:
        trimmed_contigs = Self_Align.trimmed_contigs

    }

    output{
        Array[File] trimmed_candidate_contigs = Self_Align.trimmed_contigs
    }
}


task Filter_Contigs {
    input {
        File assembly_fasta
    }

    parameter_meta {
        assembly_fasta: "Hifiasm Assembly Fasta File"
    }

    command <<<
        set -euxo pipefail
        awk 'BEGIN {RS = ">" ; ORS = ""} length($2) >16000 && length($2) < 30000 {print ">"$0}' ~{assembly_fasta} > filtered_contigs.fasta

    >>>

    output {
        File filtered_contigs = "filtered_contigs.fasta"
    }

    #########################
    runtime {
        disks: "local-disk 100 HDD"
        docker: "gcr.io/cloud-marketplace/google/ubuntu2004:latest"
    }
}


task Self_Align {
    input {
        File filtered_contigs
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        filtered_contigs: "Filtered contigs based on genome ßlength"

    }



    command <<<
        python <<CODE
        import mappy as mp
        with open("~{filtered_contigs}", "r") as f:
            assembly = f.readlines()
            for i in assembly:
                if i.startswith(">"):
                    contig_name = i.strip().split('>')[1]
                else:
                    seqlen = len(i)
                    split_pos = seqlen // 2
                    left_half = i[0:split_pos]
                    right_half = i[split_pos::]
                    left_header = contig_name+'_left'
                    right_header = contig_name+'_right'


                    with open(left_header+'.fa','w') as lfa:
                        lfa.write('>'+left_header)
                        lfa.write('\n')
                        lfa.write(left_half)
                        lfa.close()

                    with open(right_header+'.fa','w') as rfa:
                        rfa.write('>'+right_header)
                        rfa.write('\n')
                        rfa.write(right_half)
                        rfa.close()

                    aligner = mp.Aligner(left_header+'.fa')
                    if not aligner: raise Exception("ERROR: failed to load/build index")
                    for name, seq, qual in mp.fastx_read(right_header+'.fa'):
                        for hit in aligner.map(seq):
                            if hit:
                                soft_clipped = hit.q_st
                                trimmed_seq = left_half + right_half[0:hit.q_st]

                                with open(contig_name+'_trimmed.fa', 'w') as tfa:
                                    tfa.write('>'+contig_name+'_trimmed')
                                    tfa.write('\n')
                                    tfa.write(trimmed_seq)
                                    tfa.close()
        CODE
    >>>

    output {
        Array[File] trimmed_contigs = glob("*_trimmed.fa")
    }

#     Int disk_size = 2*ceil(size(filtered_contigs, "GB"))

    ########################
    #########################
#    RuntimeAttr default_attr = object {
#        cpu_cores:          1,
#        mem_gb:             8,
#        disk_gb:            disk_size,
#        boot_disk_gb:       10,
#        preemptible_tries:  1,
#        max_retries:        0,
#        docker:             "us.gcr.io/broad-dsp-lrma/lr-c3poa:2.2.2"
#    }
#    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
#    runtime {
#        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
#        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
#        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
#        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
#        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
#        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
#        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
#    }

    runtime {
        disks: "local-disk 100 HDD"
        docker: "us.gcr.io/broad-dsp-lrma/lr-c3poa:2.2.2"
    }



}




task Align_to_Candidates {

    input{
        Array[File] trimmed_contigs
        File bam
    }

    parameter_meta {

    }

    command <<<
    >>>

    output {
        Array[File] candidates_alignment_results = glob("*.bam")
    }

    runtime {
        disks: ""
        docker: ""
    }
}


