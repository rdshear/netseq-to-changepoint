version 1.0

workflow netsq_to_changepoint {

    meta {
        description: "Discover changepoints within genes in NETS-seq occupancy"
        author: "Robert D. Shear"
        email:  "rshear@gmail.com"
    }

    input {

        String sampleName
        String Cp_algorithm
        File genelist
        File CoverageBedgraph_Pos
        File CoverageBedgraph_Neg
        Int ShardCount
        Int MaxGenes
        Int GeneTrimLength
        Int MaxK 

        # environment
        Int threads = 8
        Int preemptible = 1
        String memory = "8GB"
        String docker_netcpa = 'rdshear/netcpa'
    }

    String results_filename = "~{sampleName}.cp.~{Cp_algorithm}.gff3"

    call DiscoverBreakpoints {
        input:
            genelist = genelist,
            GeneTrimLength = GeneTrimLength,
            MaxGenes = MaxGenes,
            MaxK = MaxK,
            Cp_algorithm = Cp_algorithm,
            CoverageBedgraph_Pos = CoverageBedgraph_Pos,
            CoverageBedgraph_Neg = CoverageBedgraph_Neg,
            Output_Filename = results_filename,

            docker = docker_netcpa,
            threads = threads,
            preemptible = preemptible,
            memory = memory
    }

    # call CreateShards {
    #     input:
    #         genelist = genelist,
    #         CoverageBedgraph_Pos = CoverageBedgraph_Pos,
    #         CoverageBedgraph_Neg = CoverageBedgraph_Neg,
    #         GeneTrimLength = GeneTrimLength,
    #         ShardCount = ShardCount,
    #         MaxGenes = MaxGenes,
    #         docker = docker_netcpa
    # }

    # scatter (genespec in CreateShards.shard_specs) {
    #     String Ofile = 'cp_' + basename(genespec)
    #     call DiscoverBreakpoints {
    #         input: 
    #             workfile = genespec,
    #             Output_Filename = Ofile,
    #             Cp_algorithm = Cp_algorithm,
    #             MaxK = MaxGenes,

    #             docker = docker_netcpa,
    #             threads = threads,
    #             preemptible = preemptible
    #     }
    # }

    # call GatherShards {
    #     input:
    #         outFileName = results_filename,
    #         shardResults = DiscoverBreakpoints.result_file,
    #         docker = docker_netcpa
    # }
    
    # output {
    #     File changepoint_segments = GatherShards.results
    # }
    output {
        File changepoint_segments = DiscoverBreakpoints.results
    }
}

task CreateShards {
    input {
    File genelist
    File CoverageBedgraph_Pos
    File CoverageBedgraph_Neg
    Int GeneTrimLength
    Int ShardCount
    Int MaxGenes

    String docker
    }

    command <<<
        set -e

        Rscript --vanilla /scripts/DiscoverBreakpointsScatter.R \
            ~{genelist} \
            ~{CoverageBedgraph_Pos} \
            ~{CoverageBedgraph_Neg} \
            ~{GeneTrimLength} \
            ~{MaxGenes} \
            ~{ShardCount}
    >>>

    runtime {
        docker: docker
    }

    output {
        Array[File] shard_specs = glob("shard_*.rds")
    }
}

task DiscoverBreakpoints {
    input {
        File genelist
        Int GeneTrimLength
        Int MaxGenes
        Int MaxK
        String Cp_algorithm
        File CoverageBedgraph_Pos
        File CoverageBedgraph_Neg
        String Output_Filename

        String docker
        Int threads
        Int preemptible
        String memory
    }

    command <<<
        set -e

        Rscript --vanilla /scripts/DiscoverBreakpoints2.R \
            ~{genelist} \
            ~{GeneTrimLength} \
            ~{MaxK} \
            ~{MaxGenes} \
            ~{Cp_algorithm} \
            ~{CoverageBedgraph_Pos} \
            ~{CoverageBedgraph_Neg} \
            ~{Output_Filename}
    >>>

    output {
        File results = Output_Filename
    }

    runtime {
        docker: docker
        preemptible: preemptible
        cpu: threads
        memory: memory
    }
}

task GatherShards {
    input {
        String outFileName
        Array[File] shardResults
        String docker
    }

    command <<<
        set -e

        Rscript --vanilla /scripts/DiscoverBreakpointsGather.R \
        ~{outFileName} \
        ~{sep=" " shardResults}
    >>>

    output {
        File results = outFileName
    }

    runtime {
        docker: docker
    }
}