## Overview

SampleMine is a general-purpose system for subgraph pattern mining based on subgraph enumeration and sampling. Example tasks that SampleMine supports are:

* **Subgraph Counting**: Counting the embeddings of different subgraph patterns and find the patterns with the largest counts. 
* **Frequent Subgraph Mining**: Obtaining all frequent subgraph patterns from a labeled input graph based on MNI support. 
* **User-Defined Queries**: Finding subgraphs that meet the constraints specified by the users. 


## To reproduce the results in the paper
We offered a docker image contains the executables and source code of SampleMine, Peregrine and Automine.
```shell
#The image link:
https://hub.docker.com/r/weiyihua/samplemine_test

#The code in the docker container is located in the /home directory. 
cd /home
```

The directory pattern_mining in main directory includes the instructions for testing SampleMine
The directory ComparedSystem in main directory includes two subdirectories, each of them has a README file, one for testing Automine and another for testing Peregrine. 
Please follow the instructions in 'pattern_mining' directory and 'ComparedSystem' directory for reproducing the results in the paper. 






