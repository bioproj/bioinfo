docker pull biocontainers/fastqc:v0.11.9_cv8 
docker tag biocontainers/fastqc:v0.11.9_cv8 registry.cn-hangzhou.aliyuncs.com/wybioinfo/fastqc:v0.11.9
docker push registry.cn-hangzhou.aliyuncs.com/wybioinfo/fastqc:v0.11.9


docker pull biocontainers/salmon:v0.12.0ds1-1b1-deb_cv1
docker tag biocontainers/salmon:v0.12.0ds1-1b1-deb_cv1 registry.cn-hangzhou.aliyuncs.com/wybioinfo/salmon:v0.12.0
docker push registry.cn-hangzhou.aliyuncs.com/wybioinfo/salmon:v0.12.0


docker pull nanozoo/fastp:0.23.1--9f2e255
docker tag nanozoo/fastp:0.23.1--9f2e255 registry.cn-hangzhou.aliyuncs.com/wybioinfo/fastp:0.23.1
docker push registry.cn-hangzhou.aliyuncs.com/wybioinfo/fastp:0.23.1



docker pull biocontainers/fastp:0.23.4--h5f740d0_0
docker tag biocontainers/fastp:0.23.4--h5f740d0_0 registry.cn-hangzhou.aliyuncs.com/wybioinfo/fastp:0.23.4
docker push registry.cn-hangzhou.aliyuncs.com/wybioinfo/tidyverse:4.3_1.3.0


docker pull trinityrnaseq/trinityrnaseq
docker tag trinityrnaseq/trinityrnaseq registry.cn-hangzhou.aliyuncs.com/wybioinfo/trinityrnaseq

docker pull trinityrnaseq/trinotate
docker tag trinityrnaseq/trinotate registry.cn-hangzhou.aliyuncs.com/wybioinfo/trinotate


docker pull trinityrnaseq/transdecoder
docker tag trinityrnaseq/transdecoder registry.cn-hangzhou.aliyuncs.com/wybioinfo/transdecoder
