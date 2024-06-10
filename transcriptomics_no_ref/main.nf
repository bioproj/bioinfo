#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


// params user should modify at runtime
params.reads = 'reads_{1,2}.fq.gz'


// internal parameters (typically do not need editing)
params.outprefix = './'
params.taskoutdir = 'trinity_out'
//
params.overoutdirprefix = 'overlays_'
params.overfileprefix = 'overlay_'


process overlay_one {
  tag "${dir}/${name}"
  storeDir "${dir}/${params.overoutdirprefix}${name}"

  container ''

  input:
  tuple val(dir), val(name), path(read1), path(read2)

  output:
  tuple val(dir), val(name), path("${params.overfileprefix}one")

  script:
  """
  singularity exec docker://ubuntu:18.04 bash -c ' \
  out_file=\"${params.overfileprefix}one\" && \
  mkdir -p overlay_tmp/upper overlay_tmp/work && \
  dd if=/dev/zero of=\${out_file} count=${params.overlay_size_mb_one} bs=1M && \
  mkfs.ext3 -d overlay_tmp \${out_file} && \
  rm -rf overlay_tmp \
  '
  """
}


process overlay_many {
  tag "${dir}/${name}"
  storeDir "${dir}/${params.overoutdirprefix}${name}"

  container ''

  input:
  tuple val(dir), val(name), path(reads_fa)

  output:
  tuple val(dir), val(name), val(label_fa), path("${params.overfileprefix}${label_fa}")

  script:
  label_fa = file(reads_fa).getSimpleName()
  """
  singularity exec docker://ubuntu:18.04 bash -c ' \
  out_file=\"${params.overfileprefix}${reads_fa.toString().minus('.tgz')}\" && \
  mkdir -p overlay_tmp/upper overlay_tmp/work && \
  dd if=/dev/zero of=\${out_file} count=${params.overlay_size_mb_many} bs=1M && \
  mkfs.ext3 -d overlay_tmp \${out_file} && \
  rm -rf overlay_tmp \
  '
  """
}


process jellyfish {
  // tag "${dir}/${name}"
  // stageInMode { params.copyinput ? 'copy' : 'symlink' }
  container "registry.cn-hangzhou.aliyuncs.com/wybioinfo/trinityrnaseq"
  
  input:
  tuple val(dir), val(name), path(read1), path(read2)

  output:
  tuple val(dir), val(name), path(read1), path(read2)

  script:
  """
  mem='${task.memory}'
  mem=\${mem%B}
  mem=\${mem// /}

  Trinity \
    --left $read1 \
    --right $read2 \
    --seqType fq \
    --no_normalize_reads \
    --verbose \
    --no_version_check \
    --output ${params.taskoutdir} \
    --max_memory \${mem} \
    --CPU ${task.cpus} \
    --no_run_inchworm
  """
}


process inchworm {
  // tag "${dir}/${name}"
  container "registry.cn-hangzhou.aliyuncs.com/wybioinfo/trinityrnaseq"
  input:
  tuple val(dir), val(name), path(read1), path(read2)

  output:
  tuple val(dir), val(name), path(read1), path(read2)

  script:
  """
  mem='${task.memory}'
  mem=\${mem%B}
  mem=\${mem// /}

  Trinity \
    --left $read1 \
    --right $read2 \
    --seqType fq \
    --no_normalize_reads \
    --verbose \
    --no_version_check \
    --output ${params.taskoutdir} \
    --max_memory \${mem} \
    --CPU ${task.cpus} \
    --inchworm_cpu ${task.cpus} \
    --no_run_chrysalis
  """
}


process chrysalis {
  // tag "${dir}/${name}"
  container "registry.cn-hangzhou.aliyuncs.com/wybioinfo/trinityrnaseq"

  input:
  tuple val(dir), val(name), path(read1), path(read2)

  output:
  tuple val(dir), val(name), path{"${params.taskoutdir}/read_partitions/**inity.reads.fa" }

  script:
  """


  mem='${task.memory}'
  mem=\${mem%B}
  mem=\${mem// /}

  Trinity \
    --left $read1 \
    --right $read2 \
    --seqType fq \
    --no_normalize_reads \
    --verbose \
    --no_version_check \
    --output ${params.taskoutdir} \
    --max_memory \${mem} \
    --CPU ${task.cpus} \
    --no_distributed_trinity_exec


  """
}


process butterfly {
  // tag "${dir}/${name}"
  container "registry.cn-hangzhou.aliyuncs.com/wybioinfo/trinityrnaseq"
  input:
  tuple val(dir), val(name), path(reads_fa)

  output:
  tuple val(dir), val(name), path{ "*inity.fasta" }, optional: true

// this one has been reworded compared to SIH original, and checked against Trinity code
  script:
  """
  mem='${task.memory}'
  mem=\${mem%B}
  mem=\${mem// /}


  Trinity \
    --single ${reads_fa} \
    --run_as_paired \
    --seqType fa \
    --verbose \
    --no_version_check \
    --workdir trinity_workdir \
    --output ${reads_fa}.out \
    --max_memory \${mem} \
    --CPU ${task.cpus}  \
    --trinity_complete \
    --full_cleanup \
    --no_distributed_trinity_exec

  """
}


process aggregate {
  // tag "${dir}/${name}"
  publishDir "${dir}/${params.outprefix}", mode: 'copy', saveAs: { filename -> "${name}_"+filename }
  container "registry.cn-hangzhou.aliyuncs.com/wybioinfo/trinityrnaseq"
  input:
  tuple val(dir), val(name), path(reads_fasta)

  output:
  tuple val(dir), val(name), path("Trinity.fasta"), path("Trinity.fasta.gene_trans_map")

  script:
  """
  my_trinity=\$(which Trinity)
  my_trinity=\$(dirname \$my_trinity)

  ls *inity.fasta >input_list
  cat input_list | \${my_trinity}/util/support_scripts/partitioned_trinity_aggregator.pl \
    --token_prefix TRINITY_DN --output_prefix Trinity.tmp
  mv Trinity.tmp.fasta Trinity.fasta

  \${my_trinity}/util/support_scripts/get_Trinity_gene_to_trans_map.pl Trinity.fasta > Trinity.fasta.gene_trans_map

  """
}



workflow {
    read_ch = channel.fromFilePairs( params.reads )
                .map{ it -> [ it[1][0].parent, it[0], it[1][0], it[1][1] ] }
    // read_ch.view()
    jellyfish(read_ch)
    // jellyfish.out.view()
    inchworm(jellyfish.out)

    // inchworm.out.view()
    chrysalis(inchworm.out)
    // chrysalis.out.transpose().view()
    butterfly(chrysalis.out.transpose())
    
  // butterfly.out
  //       .groupTuple(by: [0,1])
  //        .map{ yit -> [ yit[0], yit[1], yit[2].flatten() ] }.view()
    aggregate( butterfly.out
        .groupTuple(by: [0,1])
         .map{ yit -> [ yit[0], yit[1], yit[2].flatten() ] })

    // butterfly(chrysalis.out
    //   .map{ zit -> [ zit[0], zit[1], zit[2].collate( 3 ) ] }
    //   .transpose())
  // chrysalis.out
  //       .map{ zit -> [ zit[0], zit[1], zit[2] ] }
  //       .transpose()
  //       .map{ xit -> ( xit << dummy_ov) }.view()


    // chrysalis.out.view()
    // butterfly( chrysalis.out
    //   .map{ zit -> [ zit[0], zit[1], zit[2].collate( params.bf_collate ) ] }
    //   .transpose()
    //   .map{ xit -> ( xit << dummy_ov) } )

// inputs
//   if ( params.overlay ) {
//     read_ch = channel.fromFilePairs( params.reads )
//                 .map{ it -> [ it[1][0].parent, it[0], it[1][0], it[1][1] ] }
//   } else {
//     dummy_ov = file('dummy_overlay')
//     dummy_ov.text = 'dummy_overlay\n'
//     read_ch = channel.fromFilePairs( params.reads )
//                 .map{ it -> [ it[1][0].parent, it[0], it[1][0], it[1][1], dummy_ov ] }
//   }

// //
// // tasks
// //
//   if ( params.overlay ) {
//     overlay_one(read_ch)
//     jellyfish( read_ch.join(overlay_one.out, by: [0,1]) )
//   } else {
//     jellyfish(read_ch)
//   }

//   inchworm(jellyfish.out)

//   chrysalis(inchworm.out)

//   if ( params.overlay ) {
//     overlay_many( chrysalis.out.transpose() )
//     butterfly( chrysalis.out
//       .transpose().map{ xit -> [ xit[0], xit[1], file(xit[2]).getSimpleName(), xit[2] ] }
//       .join(overlay_many.out, by: [0,1,2])
//       .map{ wit -> [ wit[0], wit[1], wit[3], wit[4] ] } )
//   } else if ( params.localdisk ) {
//     butterfly( chrysalis.out
//       .transpose()
//       .map{ xit -> ( xit << dummy_ov) } )
//   } else {
//     butterfly( chrysalis.out
//       .map{ zit -> [ zit[0], zit[1], zit[2].collate( params.bf_collate ) ] }
//       .transpose()
//       .map{ xit -> ( xit << dummy_ov) } )
//   }

//   if ( params.overlay ) {
//     aggregate( butterfly.out
//       .groupTuple(by: [0,1])
//       .map{ yit -> [ yit[0], yit[1], yit[2].flatten() ] }
//       .join(overlay_one.out, by: [0,1]) )
//   } else {
//     aggregate( butterfly.out
//       .groupTuple(by: [0,1])
//       .map{ yit -> [ yit[0], yit[1], yit[2].flatten(), dummy_ov ] } )
//   }

}

// nf  run main.nf  --reads='RNASEQ_data/*.{left,right}.fq.gz' 