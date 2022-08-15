The source code of Peregrine is obtained from https://github.com/pdclab/peregrine

## Run docker image and reproduce the result of Peregrine in Table 3 and Table 6.
    
```shell
#The docker image link
https://hub.docker.com/r/weiyihua/samplemine_test

docker pull weiyihua/samplemine_test        #Download the last version of the docker image. 
docker images                               #Check image id.
docker run -it -u root xxxxxxxx /bin/bash   #Replace xxxxxxxx with the real image id. 

cd /home/peregrine

#Enable some enviroment variables
source tbb2020/bin/tbbvars.sh intel64

#Then run the following commands to do the test. 
```

```shell
#Reproduce the results of Table 3
#The meaning of the parameter is as follows:
#2nd: input graph
#3rd: pattern size
#4th minimum frequency
#5th v: vertex-induced, e: edge-induced. 
#Our subgraph counting task is same with the vertex-induced based fsm in Peregrine.
./bin/fsm data/citeseer 4 1 16 v
./bin/fsm data/citeseer 5 1 16 v
./bin/fsm data/citeseer 6 1 16 v
./bin/fsm data/citeseer 7 1 16 v

#Reproduce the results of Table 6
./bin/fsm data/citeseer 4 0.001 16 e
./bin/fsm data/citeseer 4 0.005 16 e
./bin/fsm data/citeseer 4 0.01 16 e
./bin/fsm data/citeseer 4 0.05 16 e
```
 
