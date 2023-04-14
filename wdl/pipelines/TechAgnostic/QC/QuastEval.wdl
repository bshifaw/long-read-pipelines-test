version 1.0

import "../../../tasks/QC/Quast.wdl" as Quast

workflow QuastEval {
    input {
        File reference
        Array[File] assemblies
        Array[String]? labels

        Boolean icarus = false
    }

    call Quast.QuastBenchmark as Benchmark {
        input:
            ref = reference,
            icarus = icarus,

            assemblies = assemblies,
            labels = labels
    }

    output {
        File report_html = Benchmark.report_html
        File report_txt = Benchmark.report_txt

        Array[File] plots = Benchmark.plots
        Array[File] icarus_html = Benchmark.icarus_html
    }
}
