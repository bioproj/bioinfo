#/bin/bash
docker run --rm \
    -w $PWD  \
    -v $PWD:$PWD \
    registry.cn-hangzhou.aliyuncs.com/wybioinfo/nextflowdev:23.11.0  \
    nf $@