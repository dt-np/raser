@REM # Main driver to run docker image for RASER 
@REM # Author SHI Xin <shixin@ihep.ac.cn>
@REM # Created [2021-03-31 Wed 10:56]

docker run --rm -it -h raser -e DISPLAY=192.168.237.32:0  --mount type=bind,source=C:\\Users\\dell\\raser,target=/home/physicist raser bash  







