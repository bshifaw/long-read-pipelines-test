version 1.0

import "tasks/NanoPlot.wdl"

workflow HifiUstats {

    input {
        File uBAM
    }

    call NanoPlot.NanoPlotFromUnAligned as Metrics {input: unaligned_file = uBAM, format='ubam'}

    output {
        Map[String, Float] hifi_stats_map = Metrics.stats_map
    }
}