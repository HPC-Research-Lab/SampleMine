# Peregrine

## Run docker image and reproduce the result of Peregrine in Table 3 and Table 6.
    
```shell
docker pull weiyihua/samplemine_test

#Replace 665014af3aca with the real image id. 
docker run -it -u root 665014af3aca /bin/bash

cd /home/peregrine

#Then run the following commands to do the test. 
```
```shell
#Reproduce the results of Table 3
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
 