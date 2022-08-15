## Overview

SampleMine is a general-purpose system for subgraph pattern mining based on subgraph enumeration and sampling. Example tasks that SampleMine supports are:

* **Subgraph Counting**: Counting the embeddings of different subgraph patterns and find the patterns with the largest counts. 
* **Frequent Subgraph Mining**: Obtaining all frequent subgraph patterns from a labeled input graph based on MNI support. 
* **User-Defined Queries**: Finding subgraphs that meet the constraints specified by the users. 


## To reproduce the results in the paper
We provide a docker image containing the source code and compiled executables of SampleMine, Peregrine, Automine and Pangolin.
```shell
#The image link:
https://hub.docker.com/r/weiyihua/samplemine_test


#The execution commands
docker pull weiyihua/samplemine_test        #Download the last version of the docker image. 
docker images                               #Check image id.
docker run -it -u root xxxxxxxx /bin/bash   #Replace xxxxxxxx with the real image id. 

#The code in the docker container is located in the /home directory. 
```

The 'pattern_mining' directory includes the instructions for testing SampleMine. 
The 'ComparedSystems' directory includes three subdirectories, each has a README file with the instructions for testing Automine, Peregrine and Pangolin. 





