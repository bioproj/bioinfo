params.str = 'Hello world!'

process splitLetters {
    output:
    path 'chunk_*'

    """
    printf '${params.str}' | split -b 6 - chunk_
    """
}

process convertToUpper {
    input:
    path x

    output:
    stdout

    """
    cat $x | tr '[a-z]' '[A-Z]'
    """
}
process FASTQC {
    tag "$meta.id"
    label 'process_medium'
    // https://www.nextflow.io/docs/edge/config.html#scope-docker
    // https://www.nextflow.io/docs/edge/process.html#container
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0' :
    //     'biocontainers/fastqc:0.12.1--hdfd78af_0' }"
    container "registry.cn-hangzhou.aliyuncs.com/wybioinfo/fastqc:v0.11.9"
    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip") , emit: zip
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Make list of old name and new name pairs to use for renaming in the bash while loop
    def old_new_pairs = reads instanceof Path || reads.size() == 1 ? [[ reads, "${prefix}.${reads.extension}" ]] : reads.withIndex().collect { entry, index -> [ entry, "${prefix}_${index + 1}.${entry.extension}" ] }
    def rename_to = old_new_pairs*.join(' ').join(' ')
    def renamed_files = old_new_pairs.collect{ old_name, new_name -> new_name }.join(' ')

    // def memory_in_mb = MemoryUnit.of("${task.memory}").toUnit('MB')
    // FastQC memory value allowed range (100 - 10000)
    // def fastqc_memory = memory_in_mb > 10000 ? 10000 : (memory_in_mb < 100 ? 100 : memory_in_mb)
    def fastqc_memory = 1024
    """
    printf "%s %s\\n" $rename_to | while read old_name new_name; do
        [ -f "\${new_name}" ] || ln -s \$old_name \$new_name
    done

    fastqc \\
        $args \\
        --threads $task.cpus \\
        --memory $fastqc_memory \\
        $renamed_files

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed '/FastQC v/!d; s/.*v//' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.html
    touch ${prefix}.zip

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed '/FastQC v/!d; s/.*v//' )
    END_VERSIONS
    """
}

if (!params.input) {
    exit 1, '请输入fastq路径参数,例如(/my/data/SRR*_{1,2}.fastq)!' 
} 
// input = params.input
// if(!params.input.startsWith("/")){
//     params.input="${projectDir}/${params.input}"
// }
println "##############################"
println "fastq路径参数:"+params.input
println "##############################"


workflow {
    ch_input = channel.fromFilePairs(params.input)
    // [SRR493366, [/my/data/SRR493366_1.fastq, /my/data/SRR493366_2.fastq]]
    // [[id:"SRR493366"], ["/home/wy/workspace/transcriptomics/testData/RNA-seq/reads/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz","/home/wy/workspace/transcriptomics/testData/RNA-seq/reads/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz"]]
    ch_input.map {
            meta, fastq ->
                new_id = meta - ~/$params.cutStr/
                [ [id: new_id], fastq ]
    }.set{ch_input}
    ch_input.view()
    FASTQC(ch_input)
    // splitLetters | flatten | convertToUpper | view { it.trim() }
}

