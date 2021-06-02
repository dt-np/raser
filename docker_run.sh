#!/usr/bin/env bash        

# Docker run raser image     
# Author SHI Xin <shixin@ihep.ac.cn>  
# Created [2021-03-31 Wed 10:40] 


echo "Docker run raser ..."
if [[ "$OSTYPE" == "darwin"* ]]; then 
    echo "Mac"
    xhost + 
    docker run --rm -it -h raser  -v /tmp/.X11-unix:/tmp/.X11-unix -e DISPLAY=host.docker.internal:0.0 --mount type=bind,source=$HOME/raser,target=/home/physicist raser bash 
else
        if [[ "$OSTYPE" == "linux"* ]]; then
            echo "linux"
            xhost +
            docker run --rm -it -h raser -v /tmp/.X11-unix:/tmp/.X11-unix -e DISPLAY=$DISPLAY --mount type=bind,source=$HOME/githubs/raser_test/raser,target=/home/physicist raser bash
        else
            echo "Unknow ostype:" $OSTYPE
        fi
fi  


