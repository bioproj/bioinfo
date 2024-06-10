 ${TRINITY_HOME}/Trinity --seqType fq --SS_lib_type RF  \
           --left RNASEQ_data/Sp_log.left.fq.gz,RNASEQ_data/Sp_hs.left.fq.gz,RNASEQ_data/Sp_ds.left.fq.gz,RNASEQ_data/Sp_plat.left.fq.gz \
           --right RNASEQ_data/Sp_log.right.fq.gz,RNASEQ_data/Sp_hs.right.fq.gz,RNASEQ_data/Sp_ds.right.fq.gz,RNASEQ_data/Sp_plat.right.fq.gz \
           --CPU 2 --max_memory 1G

$TRINITY_HOME/util/TrinityStats.pl trinity_out_dir/Trinity.fasta

docker run --rm  -it \
    --user $(id -u):$(id -g) \
    -v $PWD:$PWD \
    -w $PWD \
    trinityrnaseq/trinityrnaseq bash


cd /home/wy/workspace/transcriptomics/transcriptomics_no_ref/database
docker run --rm  -it \
    -v $PWD:/data \
    -v /tmp:/tmp \
    -e TRINOTATE_HOME=/usr/local/src/Trinotate \
    registry.cn-hangzhou.aliyuncs.com/wybioinfo/trinotate bash

$TRINOTATE_HOME/Trinotate  --create \
                          --db /data/myTrinotate.sqlite \
                          --trinotate_data_dir /data/db  \
                          --use_diamond | tee trinotate.log


docker run --rm  -it \
    -v $PWD:/data \
    -v /tmp:/tmp \
    registry.cn-hangzhou.aliyuncs.com/wybioinfo/transdecoder bash

TransDecoder.LongOrfs -t ../trinity_out_dir/Trinity.fasta
TransDecoder.Predict -t ../trinity_out_dir/Trinity.fasta



$TRINOTATE_HOME/Trinotate --db myTrinotate.sqlite --init \
           --gene_trans_map ../trinity_out_dir/Trinity.fasta.gene_trans_map \
           --transcript_fasta ../trinity_out_dir/Trinity.fasta \
           --transdecoder_pep ../transdecoder/Trinity.fasta.transdecoder.pep

$TRINOTATE_HOME/Trinotate --db myTrinotate.sqlite --CPU 3 \
           --transcript_fasta ../trinity_out_dir/Trinity.fasta \
           --transdecoder_pep ../transdecoder/Trinity.fasta.transdecoder.pep \
            --trinotate_data_dir ../database/db \
            --run "swissprot_blastp swissprot_blastx pfam signalp6 tmhmmv2 infernal EggnogMapper" \
            --use_diamond 


$TRINOTATE_HOME/Trinotate --db myTrinotate.sqlite --incl_pep --report > myTrinotate.incl_pep.tsv


$TRINOTATE_HOME/util/Trinotate_get_feature_name_encoding_attributes.pl \
                  myTrinotate.tsv  > annot_feature_map.txt