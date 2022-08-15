# Automine

## Run docker image and reproduce the result of Automine in Table 3 and Table 6.
    
```shell
docker pull weiyihua/samplemine_test

#Replace 665014af3aca with the real image id. 
docker run -it -u root 665014af3aca /bin/bash

cd /home/SampleMine/pattern_mining

#Then run the following commands to do the test. 
```
```shell
#Reproduce the results of Table 3
./automine_sc.exe ./data/citeseer.lg 4
./automine_sc.exe ./data/citeseer.lg 5
./automine_sc.exe ./data/citeseer.lg 6
./automine_sc.exe ./data/citeseer.lg 7

#Reproduce the results of Table 6
./automine_fsm.exe ./data/citeseer.lg 4
./automine_fsm.exe ./data/citeseer.lg 5
./automine_fsm.exe ./data/citeseer.lg 6
./automine_fsm.exe ./data/mico.lg 4
```
 