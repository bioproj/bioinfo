#/bin/bash
docker run --rm  -it \
    -v /var/run/docker.sock:/var/run/docker.sock \
    -v  /usr/bin/docker:/usr/bin/docker \
    -v /etc/docker/daemon.json:/etc/docker/daemon.json \
    -v /var/lib/docker:/var/lib/docker \
    -w $PWD  \
    -v $PWD:$PWD  \
    -v /home/wy/workspace/transcriptomics:/home/wy/workspace/transcriptomics \
    -v /home/wy/workspace/genomics:/home/wy/workspace/genomics \
    -e NXF_LOG_FILE=logs/.nextflow.log \
    -e NXF_USRMAP=$(id -u) \
    -e NXF_GROUPMAP=$(getent group docker | cut -d: -f3) \
    registry.cn-hangzhou.aliyuncs.com/wybioinfo/nextflowdev:23.11.0  \
    nf $@ 