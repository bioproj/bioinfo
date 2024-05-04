[ ! -d  testData ] && mkdir testData
ln -s ../../testData/RNA-seq/reads/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz  testData/HBR_Rep1.read1.fastq.gz
ln -s ../../testData/RNA-seq/reads/HBR_Rep1_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz  testData/HBR_Rep1.read2.fastq.gz
ln -s ../../testData/RNA-seq/reads/HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz  testData/HBR_Rep2.read1.fastq.gz
ln -s ../../testData/RNA-seq/reads/HBR_Rep2_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz  testData/HBR_Rep2.read2.fastq.gz
ln -s ../../testData/RNA-seq/reads/HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22.read1.fastq.gz  testData/HBR_Rep3.read1.fastq.gz
ln -s ../../testData/RNA-seq/reads/HBR_Rep3_ERCC-Mix2_Build37-ErccTranscripts-chr22.read2.fastq.gz  testData/HBR_Rep3.read2.fastq.gz

ln -s ../../testData/RNA-seq/reads/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz  testData/UHR_Rep1.read1.fastq.gz
ln -s ../../testData/RNA-seq/reads/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz  testData/UHR_Rep1.read2.fastq.gz
ln -s ../../testData/RNA-seq/reads/UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz  testData/UHR_Rep2.read1.fastq.gz
ln -s ../../testData/RNA-seq/reads/UHR_Rep2_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz  testData/UHR_Rep2.read2.fastq.gz
ln -s ../../testData/RNA-seq/reads/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz  testData/UHR_Rep3.read1.fastq.gz
ln -s ../../testData/RNA-seq/reads/UHR_Rep3_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz  testData/UHR_Rep3.read2.fastq.gz


# ln -s ../../testData/RNA-seq/genomic/chr22_with_ERCC92.fa.gz  testData/chr22_with_ERCC92.fa.gz
ln -s ../../testData/RNA-seq/genomic/chr22_with_ERCC92.gtf.gz  testData/chr22_with_ERCC92.gtf.gz
ln -s ../../testData/RNA-seq/genomic/chr22_transcripts.fa.gz  testData/chr22_transcripts.fa.gz
ln -s ../../testData/RNA-seq/genomic/chr22.fa.gz  testData/chr22.fa.gz






# gzip -dc chr22_with_ERCC92.fa.gz > chr22_with_ERCC92.fa