docker pull biocontainers/fastqc:v0.11.9_cv8 
docker tag biocontainers/fastqc:v0.11.9_cv8 registry.cn-hangzhou.aliyuncs.com/wybioinfo/fastqc:v0.11.9
docker push registry.cn-hangzhou.aliyuncs.com/wybioinfo/fastqc:v0.11.9


docker pull biocontainers/salmon:v0.12.0ds1-1b1-deb_cv1
docker tag biocontainers/salmon:v0.12.0ds1-1b1-deb_cv1 registry.cn-hangzhou.aliyuncs.com/wybioinfo/salmon:v0.12.0
docker push registry.cn-hangzhou.aliyuncs.com/wybioinfo/salmon:v0.12.0


docker pull rocker/tidyverse
docker tag rocker/tidyverse registry.cn-hangzhou.aliyuncs.com/wybioinfo/tidyverse:4.3_1.3.0
docker push registry.cn-hangzhou.aliyuncs.com/wybioinfo/tidyverse:4.3_1.3.0

registry.cn-hangzhou.aliyuncs.com/wybioinfo/sra-tools:3.1.0
