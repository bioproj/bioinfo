#!/bin/bash
docker run --rm  \
    -e USERID=$(id -u) -e GROUPID=$(id -g) \
    -e DISABLE_AUTH=true \
    -e R_LIBS=$PWD/R_LIBS  \
    -w $PWD  \
    -v $PWD:$PWD  \
    -v $PWD:/home/rstudio  \
    -p 8787:8787  \
    registry.cn-hangzhou.aliyuncs.com/wybioinfo/tidyverse:4.4
