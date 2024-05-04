docker run --rm  \
    -e USERID=$(id -u) -e GROUPID=$(id -g) \
    -e DISABLE_AUTH=true \
    -e R_LIBS=/home/wy/workspace/transcriptomics/transcriptomics/R_LIBS  \
    -w $PWD  \
    -v $PWD:$PWD  \
    -v $PWD:/home/rstudio  \
    -v /home/wy/workspace/transcriptomics:/home/wy/workspace/transcriptomics \
    -v /home/wy/workspace/genomics:/home/wy/workspace/genomics \
    -p 8787:8787  \
    registry.cn-hangzhou.aliyuncs.com/wybioinfo/tidyverse:4.4
