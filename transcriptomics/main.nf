params.str = 'Hello world!'



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
process FASTP {
    tag "$meta.id"
    label 'process_medium'

    container "registry.cn-hangzhou.aliyuncs.com/wybioinfo/fastp:0.23.1"

    input:
    tuple val(meta), path(reads)
    val  adapter_fasta
    val   save_trimmed_fail
    val   save_merged

    output:
    tuple val(meta), path('*.fastp.fastq.gz') , optional:true, emit: reads
    tuple val(meta), path('*.json')           , emit: json
    tuple val(meta), path('*.html')           , emit: html
    tuple val(meta), path('*.log')            , emit: log
    path "versions.yml"                       , emit: versions
    tuple val(meta), path('*.fail.fastq.gz')  , optional:true, emit: reads_fail
    tuple val(meta), path('*.merged.fastq.gz'), optional:true, emit: reads_merged

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def adapter_list = adapter_fasta ? "--adapter_fasta ${adapter_fasta}" : ""
    def fail_fastq = save_trimmed_fail && meta.single_end ? "--failed_out ${prefix}.fail.fastq.gz" : save_trimmed_fail && !meta.single_end ? "--failed_out ${prefix}.paired.fail.fastq.gz --unpaired1 ${prefix}_1.fail.fastq.gz --unpaired2 ${prefix}_2.fail.fastq.gz" : ''
    // Added soft-links to original fastqs for consistent naming in MultiQC
    // Use single ended for interleaved. Add --interleaved_in in config.
    if (meta.single_end) {
        """
        [ ! -f  ${prefix}.fastq.gz ] && ln -sf $reads ${prefix}.fastq.gz

        fastp \\
            --in1 ${prefix}.fastq.gz \\
            --out1  ${prefix}.fastp.fastq.gz \\
            --thread $task.cpus \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html \\
            $adapter_list \\
            $fail_fastq \\
            $args \\
            2> >(tee ${prefix}.fastp.log >&2)

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
        END_VERSIONS
        """
    } else {
        def merge_fastq = save_merged ? "-m --merged_out ${prefix}.merged.fastq.gz" : ''
        """
        [ ! -f  ${prefix}_1.fastq.gz ] && ln -sf ${reads[0]} ${prefix}_1.fastq.gz
        [ ! -f  ${prefix}_2.fastq.gz ] && ln -sf ${reads[1]} ${prefix}_2.fastq.gz
        fastp \\
            --in1 ${prefix}_1.fastq.gz \\
            --in2 ${prefix}_2.fastq.gz \\
            --out1 ${prefix}_1.fastp.fastq.gz \\
            --out2 ${prefix}_2.fastp.fastq.gz \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html \\
            $adapter_list \\
            $fail_fastq \\
            $merge_fastq \\
            --thread $task.cpus \\
            --detect_adapter_for_pe \\
            $args \\
            2> >(tee ${prefix}.fastp.log >&2)

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
        END_VERSIONS
        """
    }

    stub:
    def prefix              = task.ext.prefix ?: "${meta.id}"
    def is_single_output    = task.ext.args?.contains('--interleaved_in') || meta.single_end
    def touch_reads         = is_single_output ? "${prefix}.fastp.fastq.gz" : "${prefix}_1.fastp.fastq.gz ${prefix}_2.fastp.fastq.gz"
    def touch_merged        = (!is_single_output && save_merged) ? "touch ${prefix}.merged.fastq.gz" : ""
    """
    touch $touch_reads
    touch "${prefix}.fastp.json"
    touch "${prefix}.fastp.html"
    touch "${prefix}.fastp.log"
    $touch_merged

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
    END_VERSIONS
    """
}
process SALMON_INDEX {
    tag "$transcript_fasta"
    label "process_medium"

    container "registry.cn-hangzhou.aliyuncs.com/wybioinfo/salmon:v1.10.3"

    input:
    val genome_fasta 
    path transcript_fasta 

    output:
    path "salmon"      , emit: index
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def get_decoy_ids = "grep '^>' $genome_fasta | cut -d ' ' -f 1 | cut -d \$'\\t' -f 1 > decoys.txt"

    if(!genome_fasta){
        """
        salmon index \\
            --threads $task.cpus \\
            -t $transcript_fasta \\
            -i salmon
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            salmon: \$(echo \$(salmon --version) | sed -e "s/salmon //g")
        END_VERSIONS
        """
    }else{
        def gentrome      = "gentrome.fa"
        if (genome_fasta.toString().endsWith('.gz')) {
        get_decoy_ids = "grep '^>' <(gunzip -c $genome_fasta) | cut -d ' ' -f 1 | cut -d \$'\\t' -f 1 > decoys.txt"
        gentrome      = "gentrome.fa.gz"
        }
        """
        $get_decoy_ids
        sed -i.bak -e 's/>//g' decoys.txt
        cat $transcript_fasta $genome_fasta > $gentrome

        salmon \\
            index \\
            --threads $task.cpus \\
            -t $gentrome \\
            -d decoys.txt \\
            $args \\
            -i salmon

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            salmon: \$(echo \$(salmon --version) | sed -e "s/salmon //g")
        END_VERSIONS
        """
    }


    stub:
    """
    mkdir salmon
    touch salmon/complete_ref_lens.bin
    touch salmon/ctable.bin
    touch salmon/ctg_offsets.bin
    touch salmon/duplicate_clusters.tsv
    touch salmon/info.json
    touch salmon/mphf.bin
    touch salmon/pos.bin
    touch salmon/pre_indexing.log
    touch salmon/rank.bin
    touch salmon/refAccumLengths.bin
    touch salmon/ref_indexing.log
    touch salmon/reflengths.bin
    touch salmon/refseq.bin
    touch salmon/seq.bin
    touch salmon/versionInfo.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(echo \$(salmon --version) | sed -e "s/salmon //g")
    END_VERSIONS
    """
}

process SALMON_QUANT {
    publishDir = [
        path: { "${params.outdir}/salmon" },
        mode: 'copy',
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]
    tag "$meta.id"
    label "process_medium"

    container "registry.cn-hangzhou.aliyuncs.com/wybioinfo/salmon:v1.10.3"

    input:
    tuple val(meta), path(reads)
    path  index
    path  gtf
    path  transcript_fasta
    val   alignment_mode
    val   lib_type

    output:
    tuple val(meta), path("${prefix}") , emit: results
    tuple val(meta), path("*info.json"), emit: json_info, optional: true
    path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"

    def reference   = "--index $index"
    def reads1 = [], reads2 = []
    meta.single_end ? [reads].flatten().each{reads1 << it} : reads.eachWithIndex{ v, ix -> ( ix & 1 ? reads2 : reads1) << v }
    def input_reads = meta.single_end ? "-r ${reads1.join(" ")}" : "-1 ${reads1.join(" ")} -2 ${reads2.join(" ")}"
    if (alignment_mode) {
        reference   = "-t $transcript_fasta"
        input_reads = "-a $reads"
    }

    def strandedness_opts = [
        'A', 'U', 'SF', 'SR',
        'IS', 'IU' , 'ISF', 'ISR',
        'OS', 'OU' , 'OSF', 'OSR',
        'MS', 'MU' , 'MSF', 'MSR'
    ]
    def strandedness =  'A'
    if (lib_type) {
        if (strandedness_opts.contains(lib_type)) {
            strandedness = lib_type
        } else {
            log.info "[Salmon Quant] Invalid library type specified '--libType=${lib_type}', defaulting to auto-detection with '--libType=A'."
        }
    } else {
        strandedness = meta.single_end ? 'U' : 'IU'
        if (meta.strandedness == 'forward') {
            strandedness = meta.single_end ? 'SF' : 'ISF'
        } else if (meta.strandedness == 'reverse') {
            strandedness = meta.single_end ? 'SR' : 'ISR'
        }
    }
            // --geneMap $gtf \\
    """
    salmon quant \\
        --threads $task.cpus \\
        --libType=$strandedness \\
        $reference \\
        $input_reads \\
        $args \\
        -o $prefix

    if [ -f $prefix/aux_info/meta_info.json ]; then
        cp $prefix/aux_info/meta_info.json "${prefix}_meta_info.json"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(echo \$(salmon --version) | sed -e "s/salmon //g")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    touch ${prefix}_meta_info.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(echo \$(salmon --version) | sed -e "s/salmon //g")
    END_VERSIONS
    """
}

process SALMON_COUNT_R {
    publishDir = [
        path: { "${params.outdir}/R" },
        mode: 'copy',
    ]
    debug true
    container "registry.cn-hangzhou.aliyuncs.com/wybioinfo/rstudio_transcriptomics:4.4"

    input:
    path quant_sf_dir
    path gtf
    val treatment_group
    val control_group

    output:
    path "count.tsv"

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # echo \$(which count.R )
    count.R  $gtf $quant_sf_dir $treatment_group $control_group
    """
    //    count.R  $gtf $quant_sf_dir $treatment_group $control_group

    output:
    path  "count.tsv"
    script:
    if(treatment_group && control_group){
    """
    count.R  $gtf ./ $treatment_group $control_group
    """
    }else{
    """
    # echo \$(which count.R )
    count.R  $gtf ./
    """
    }

    //    count.R  $gtf ./
    // Rscript bin/count.R  ./testData/chr22_genes.gtf out/salmon
}




if (!params.input) {
    exit 1, '请输入fastq路径参数,例如(/my/data/SRR*_{1,2}.fastq)!' 
} 

if (!params.gtf_file) {
    exit 1, '请输入gtf路径参数,例如(testData/chr22_genes.gtf)!' 
} 
gtf_file = file(params.gtf_file)
if(!params.gtf_file.startsWith("/")){
    gtf_file =file("${projectDir}/${params.gtf_file}")
}

if (!params.transcripts_fa) {
    exit 1, '请输入transcripts路径参数,例如(testData/chr22_transcripts.fa.gz)!' 
} 
transcripts_fa = file(params.transcripts_fa)
if(!params.transcripts_fa.startsWith("/")){
    transcripts_fa =file("${projectDir}/${params.transcripts_fa}")
}

genome_fa = params.genome_fa
if(!genome_fa){
    genome_fa = false
}else{
    genome_fa = file(genome_fa)
    if(!genome_fa.startsWith("/")){
        genome_fa =file("${projectDir}/${genome_fa}")
    }
}

control_group = params.control_group
if(!control_group){
    control_group = false
}
treatment_group = params.treatment_group
if(!treatment_group){
    treatment_group = false
}


println genome_fa
println "##############################"
println "fastq路径参数:"+params.input
println "gtf路径参数:"+gtf_file
println "transcripts路径参数:"+transcripts_fa
if(genome_fa){
    println "genome路径参数:"+genome_fa
}
println "##############################"


workflow {
    ch_input = channel.fromFilePairs(params.input)
    // [SRR493366, [/my/data/SRR493366_1.fastq, /my/data/SRR493366_2.fastq]]
    // [[id:"SRR493366"], ["/home/wy/workspace/transcriptomics/testData/RNA-seq/reads/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz","/home/wy/workspace/transcriptomics/testData/RNA-seq/reads/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz"]]
    // ch_input.view()
    ch_input.map {
            meta, fastq ->
                new_id = meta - ~/$params.cutStr/
                [ [id: new_id,single_end:false], fastq ]
    }.set{ch_input}
    FASTQC(ch_input)
    ch_fastp_out = FASTP(ch_input,false,false,false)
    ch_input = ch_fastp_out.reads

    ch_salman_index = SALMON_INDEX(genome_fa,transcripts_fa).index
    // ch_salman_index.view()

    ch_salmon_quant_out = SALMON_QUANT(ch_input, ch_salman_index,gtf_file,transcripts_fa,false,"A")
    ch_salmon_quant_out_collect = ch_salmon_quant_out.results.collect{it[1]}
    // ch_salmon_quant_out_collect.view()
    SALMON_COUNT_R(ch_salmon_quant_out_collect,gtf_file,treatment_group,control_group)
    // SALMON_INDEX("/home/wy/workspace/genomics/GRCm38.primary_assembly.genome.fa.gz","/home/wy/workspace/genomics/gencode.vM23.transcripts.fa.gz")
    
    // splitLetters | flatten | convertToUpper | view { it.trim() }
}

