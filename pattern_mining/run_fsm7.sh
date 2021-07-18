#!/bin/bash

./fsm7_two_vertex.exe /data/not_backed_up/shared/gpm_data/citeseer.lg 0.001 0 > output/fsm7_citeseer_001.txt
./fsm7_two_vertex.exe /data/not_backed_up/shared/gpm_data/citeseer.lg 0.005 0 > output/fsm7_citeseer_005.txt
./fsm7_two_vertex.exe /data/not_backed_up/shared/gpm_data/citeseer.lg 0.01 0 > output/fsm7_citeseer_01.txt
./fsm7_two_vertex.exe /data/not_backed_up/shared/gpm_data/citeseer.lg 0.05 0 > output/fsm7_citeseer_05.txt
