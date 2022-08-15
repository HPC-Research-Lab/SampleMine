## Building SampleMine

### Run docker image for testing SampleMine
We provide a docker image for testing SampleMine, Automine and Peregrine. 
You can download the docker image from the following link and omit the 'Install Dependencies' and 'Compile SampleMine' process.

```Shell
#The docker image link
https://hub.docker.com/r/weiyihua/samplemine_test

docker pull weiyihua/samplemine_test        #Download the last version of the docker image. 
docker images                               #Check image id.
docker run -it -u root xxxxxxxx /bin/bash   #Replace xxxxxxxx with the real image id. 

cd /home/SampleMine/pattern_mining/
```

The source code of SampleMine, Automine and Peregrine in the container are in the /home/ directory. 
We integrated the source code of Automine into SampleMine so you will see 2 directories in /home/ directory. 
The detailed instructions for testing peregrine and Automine are in READMEs under /SampleMine/ComparedSystem.

### Install Dependencies

Download Boost Library version 1.76
```Shell
wget https://boostorg.jfrog.io/artifactory/main/release/1.76.0/source/boost_1_76_0.tar.gz
tar -xzvf boost_1_76_0.tar.gz
```

### Compile SampleMine

Download Source Code
```Shell
git clone SampleMine
```
Add enviroment variable BOOST_ROOT for compiling
```Shell
export BOOST_ROOT=PATH_TO_BOOST/boost_1_76_0/
```
Enter into the directory and compile SampleMine
```Shell
cd SampleMine/pattern_mining/
make -j 16
```


## Reproducing results of Table3, Figure9, Table4, Table5, Table6, Figure10, Table 7, Figure 11, Figure 12 and Table 8.

Before start testing, export the following 2 enviroment variables first.
```Shell
export DB_PATH=./db
export OMP_NUM_THREADS=16
```

We provide two test graphs(mico and citeseer) under ./data directory in this github repo. 
The other input graphs can be downloaded from this google drive link:[Graphs](https://drive.google.com/drive/folders/1sK-JFJim1e3N1Qd5XIupXbJTHjzasCpn?usp=sharing). After downloading, you need to store these graphs in ./data directory. 

### Table3 & Figure9

```Shell
#For TV-smpl column
./sc4_tv.exe ./data/citeseer.lg 1 2 4
./sc5_tv.exe ./data/citeseer.lg 1 2 4
./sc6_tv.exe ./data/citeseer.lg 1 2 4
./sc7_tv.exe ./data/citeseer.lg 1 2 4
./sc4_tv.exe ./data/mico.lg 1 8 64

#For SV-smpl column
./sc4_sv.exe ./data/citeseer.lg 2
./sc5_sv.exe ./data/citeseer.lg 2
./sc6_sv.exe ./data/citeseer.lg 2
./sc7_sv.exe ./data/citeseer.lg 2
./sc4_sv.exe ./data/mico.lg 8

#For TV-acc column
./sc4_tv.exe ./data/citeseer.lg 1 1 1
./sc5_tv.exe ./data/citeseer.lg 1 1 1
./sc6_tv.exe ./data/citeseer.lg 1 1 1
./sc7_tv.exe ./data/citeseer.lg 1 1 1
./sc4_tv.exe ./data/mico.lg 1 1 1

#For SV-acc column
./sc4_sv.exe ./data/citeseer.lg 1
./sc5_sv.exe ./data/citeseer.lg 1
./sc6_sv.exe ./data/citeseer.lg 1
./sc7_sv.exe ./data/citeseer.lg 1
./sc4_sv.exe ./data/mico.lg 1
```

### Table 4
```Shell
#The Tot#, SM-tot and SM-50 columns can be verified by the TV-acc execution of Table3
#From each of the following execution, you can get 3 results corresponding to the Tot#  ASAP-tot ASAP-50 respectively. 
./sc4_tv_asap.exe ./data/citeseer.lg 1 2 4
./sc5_tv_asap.exe ./data/citeseer.lg 1 2 4
./sc6_tv_asap.exe ./data/citeseer.lg 1 2 4
./sc7_tv_asap.exe ./data/citeseer.lg 1 2 4
./sc4_tv_asap.exe ./data/mico.lg 1 8 64
```

### Table 5
```Shell
./sc5_tv.exe ./data/mico.lg 1 8 64
./sc5_tv.exe ./data/com-orkut.lg 8 32 1024
./sc4_tv.exe ./data/uk-2005.lg 1024 512 1024
./sc4_tv.exe ./data/com-Friendster.lg 32 8 32
```

### Table 6
```Shell
#For SM column
./fsm4_tv.exe ./data/citeseer.lg 0.001 4
./fsm4_tv.exe ./data/citeseer.lg 0.005 4
./fsm4_tv.exe ./data/citeseer.lg 0.01 4
./fsm4_tv.exe ./data/citeseer.lg 0.05 4

./fsm4_tv.exe ./data/mico.lg 0.001 6
./fsm4_tv.exe ./data/mico.lg 0.005 6
./fsm4_tv.exe ./data/mico.lg 0.01 6
./fsm4_tv.exe ./data/mico.lg 0.05 6

./fsm5_tv.exe ./data/citeseer.lg 0.001 4
./fsm5_tv.exe ./data/citeseer.lg 0.005 4
./fsm5_tv.exe ./data/citeseer.lg 0.01 4
./fsm5_tv.exe ./data/citeseer.lg 0.05 4

./fsm6_tv.exe ./data/citeseer.lg 0.001 4
./fsm6_tv.exe ./data/citeseer.lg 0.005 4
./fsm6_tv.exe ./data/citeseer.lg 0.01 4
./fsm6_tv.exe ./data/citeseer.lg 0.05 4

#For TV-acc column
./fsm4_tv.exe ./data/citeseer.lg 0.001 0
./fsm4_tv.exe ./data/citeseer.lg 0.005 0
./fsm4_tv.exe ./data/citeseer.lg 0.01 0
./fsm4_tv.exe ./data/citeseer.lg 0.05 0

./fsm4_tv.exe ./data/mico.lg 0.001 0
./fsm4_tv.exe ./data/mico.lg 0.005 0
./fsm4_tv.exe ./data/mico.lg 0.01 0
./fsm4_tv.exe ./data/mico.lg 0.05 0

./fsm5_tv.exe ./data/citeseer.lg 0.001 0
./fsm5_tv.exe ./data/citeseer.lg 0.005 0
./fsm5_tv.exe ./data/citeseer.lg 0.01 0
./fsm5_tv.exe ./data/citeseer.lg 0.05 0

./fsm6_tv.exe ./data/citeseer.lg 0.001 0
./fsm6_tv.exe ./data/citeseer.lg 0.005 0
./fsm6_tv.exe ./data/citeseer.lg 0.01 0
./fsm6_tv.exe ./data/citeseer.lg 0.05 0


#For SV-acc column
./fsm4_sv.exe ./data/citeseer.lg 0.001 0
./fsm4_sv.exe ./data/citeseer.lg 0.005 0
./fsm4_sv.exe ./data/citeseer.lg 0.01 0
./fsm4_sv.exe ./data/citeseer.lg 0.05 0

./fsm4_sv.exe ./data/mico.lg 0.001 0
./fsm4_sv.exe ./data/mico.lg 0.005 0
./fsm4_sv.exe ./data/mico.lg 0.01 0
./fsm4_sv.exe ./data/mico.lg 0.05 0

./fsm5_sv.exe ./data/citeseer.lg 0.001 0
./fsm5_sv.exe ./data/citeseer.lg 0.005 0
./fsm5_sv.exe ./data/citeseer.lg 0.01 0
./fsm5_sv.exe ./data/citeseer.lg 0.05 0

./fsm6_sv.exe ./data/citeseer.lg 0.001 0
./fsm6_sv.exe ./data/citeseer.lg 0.005 0
./fsm6_sv.exe ./data/citeseer.lg 0.01 0
./fsm6_sv.exe ./data/citeseer.lg 0.05 0
```

### Figure 10
```Shell

./fsm4_tv.exe ./data/mico.lg 0.001 0
./fsm4_tv.exe ./data/mico.lg 0.001 4
./fsm4_tv.exe ./data/mico.lg 0.001 6
./fsm4_tv.exe ./data/mico.lg 0.001 8
./fsm4_tv.exe ./data/mico.lg 0.001 10

./fsm4_tv.exe ./data/mico.lg 0.005 0
./fsm4_tv.exe ./data/mico.lg 0.005 4
./fsm4_tv.exe ./data/mico.lg 0.005 6
./fsm4_tv.exe ./data/mico.lg 0.005 8
./fsm4_tv.exe ./data/mico.lg 0.005 10

./fsm4_tv.exe ./data/mico.lg 0.01 0
./fsm4_tv.exe ./data/mico.lg 0.01 4
./fsm4_tv.exe ./data/mico.lg 0.01 6
./fsm4_tv.exe ./data/mico.lg 0.01 8
./fsm4_tv.exe ./data/mico.lg 0.01 10

./fsm4_tv.exe ./data/mico.lg 0.05 0
./fsm4_tv.exe ./data/mico.lg 0.05 4
./fsm4_tv.exe ./data/mico.lg 0.05 6
./fsm4_tv.exe ./data/mico.lg 0.05 8
./fsm4_tv.exe ./data/mico.lg 0.05 10
```

### Table 7
```Shell
#The third argument is Sup
./table7.exe ./data/uk-2005.lg 0.0001 4 2
./table7.exe ./data/uk-2005.lg 0.0001 16 2
./table7.exe ./data/uk-2005.lg 0.0001 4 4
./table7.exe ./data/uk-2005.lg 0.0001 16 4

./table7.exe ./data/uk-2005.lg 0.0005 4 2
./table7.exe ./data/uk-2005.lg 0.0005 16 2
./table7.exe ./data/uk-2005.lg 0.0005 4 4
./table7.exe ./data/uk-2005.lg 0.0005 16 4
```


### Figure 11
```Shell
#For size-7 on CI
./q1_size7.exe ./data/citeseer.lg 0
./q1_size7.exe ./data/citeseer.lg 2
./q1_size7.exe ./data/citeseer.lg 4

./q2_size7.exe ./data/citeseer.lg 0
./q2_size7.exe ./data/citeseer.lg 2
./q2_size7.exe ./data/citeseer.lg 4

./q3_size7.exe ./data/citeseer.lg 0
./q3_size7.exe ./data/citeseer.lg 2
./q4_size7.exe ./data/citeseer.lg 4

./q4_size7.exe ./data/citeseer.lg 0
./q4_size7.exe ./data/citeseer.lg 2
./q4_size7.exe ./data/citeseer.lg 4

#For size-4 on MI
./q1_size4.exe ./data/mico.lg 0
./q1_size4.exe ./data/mico.lg 32
./q1_size4.exe ./data/mico.lg 64

./q2_size4.exe ./data/mico.lg 0
./q2_size4.exe ./data/mico.lg 32
./q2_size4.exe ./data/mico.lg 64

./q3_size4.exe ./data/mico.lg 0
./q3_size4.exe ./data/mico.lg 32
./q3_size4.exe ./data/mico.lg 64

./q4_size4.exe ./data/mico.lg 0
./q4_size4.exe ./data/mico.lg 32
./q4_size4.exe ./data/mico.lg 64

```

### Figure 12
```Shell
./q5_size7.exe ./data/citeseer.lg 0.001 0
./q5_size7.exe ./data/citeseer.lg 0.001 4
./q5_size7.exe ./data/citeseer.lg 0.001 6
./q5_size7.exe ./data/citeseer.lg 0.001 8
./q5_size7.exe ./data/citeseer.lg 0.001 10

./q5_size4.exe ./data/mico.lg 0.001 0
./q5_size4.exe ./data/mico.lg 0.001 4
./q5_size4.exe ./data/mico.lg 0.001 6
./q5_size4.exe ./data/mico.lg 0.001 8
./q5_size4.exe ./data/mico.lg 0.001 10
```

### Table 8
```Shell
./table8_q1_size5_MI.exe ./data/mico.lg 64
./table8_q1_size5_OK.exe ./data/com-orkut.lg 1024
./table8_q1_size4_UK.exe ./data/uk-2005.lg 1024
./table8_q1_size4_FR.exe ./data/com-Friendster.lg 32

./table8_q2_size5_OK.exe ./data/com-orkut.lg 1024
./table8_q2_size4_UK.exe ./data/uk-2005.lg 1024
./table8_q2_size4_FR.exe ./data/com-Friendster.lg 32
```
