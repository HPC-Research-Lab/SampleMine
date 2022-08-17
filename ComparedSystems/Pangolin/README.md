The source code of Pangolin is obtained from https://github.com/IntelligentSoftwareSystems/Galois/tree/master

## Run docker image and reproduce the result of Pangolin in Table 6.
    
```shell
#The docker image link
https://hub.docker.com/r/weiyihua/samplemine_test

docker pull weiyihua/samplemine_test        #Download the last version of the docker image. 
docker images                               #Check image id.
docker run -it -u root xxxxxxxx /bin/bash   #Replace xxxxxxxx with the real image id. 

cd /home/Galois/build/lonestar/mining/cpu/frequent-subgraph-mining/

#Then run the following commands to do the test. 
```

```shell
#Reproduce the results of Table 6
#The meaning of the parameters is as follows:
#2nd: related to the graph attribute
#3rd: related to the graph attribute
#4th: the input graph
#5th: graph format option
#6th: graph format is `adj`
#7th: pattern size: 4.
#8th: minimum support, corresponding to the for rows in the table 6
#9th: -t 1: number of trials
#Our subgraph counting task is same with the vertex-induced based fsm in Pangolin.
./frequent-subgraph-mining-cpu -symmetricGraph -simpleGraph citeseer.adj -ft adj -k=4 -ms=1 -t 1
./frequent-subgraph-mining-cpu -symmetricGraph -simpleGraph citeseer.adj -ft adj -k=4 -ms=5 -t 1
./frequent-subgraph-mining-cpu -symmetricGraph -simpleGraph citeseer.adj -ft adj -k=4 -ms=10 -t 1
./frequent-subgraph-mining-cpu -symmetricGraph -simpleGraph citeseer.adj -ft adj -k=4 -ms=50 -t 1


```
 
