#!/bin/bash

echo "Started" | mail -s "Job" peng-jiang@uiowa.edu


./fsm_ns.exe /data/not_backed_up/shared/graphs/citeseer.lg 0.001 1> fsm_ns_citeseer_7_001.output 2> fsm_ns_citeseer_7_001.err 

./fsm_ns.exe /data/not_backed_up/shared/graphs/citeseer.lg 0.005 1> fsm_ns_citeseer_7_005.output 2> fsm_ns_citeseer_7_005.err

./fsm_ns.exe /data/not_backed_up/shared/graphs/citeseer.lg 0.01 1> fsm_ns_citeseer_7_01.output 2> fsm_ns_citeseer_7_01.err

./fsm_ns.exe /data/not_backed_up/shared/graphs/citeseer.lg 0.05 1> fsm_ns_citeseer_7_05.output 2> fsm_ns_citeseer_7_05.err

echo "Done" | mail -s "Job1" peng-jiang@uiowa.edu

./fsm8_ns.exe /data/not_backed_up/shared/graphs/citeseer.lg 0.001 1> fsm_ns_citeseer_8_001.output 2> fsm_ns_citeseer_8_001.err 

./fsm8_ns.exe /data/not_backed_up/shared/graphs/citeseer.lg 0.005 1> fsm_ns_citeseer_8_005.output 2> fsm_ns_citeseer_8_005.err

./fsm8_ns.exe /data/not_backed_up/shared/graphs/citeseer.lg 0.01 1> fsm_ns_citeseer_8_01.output 2> fsm_ns_citeseer_8_01.err

./fsm8_ns.exe /data/not_backed_up/shared/graphs/citeseer.lg 0.05 1> fsm_ns_citeseer_8_05.output 2> fsm_ns_citeseer_8_05.err

echo "Done" | mail -s "Job2" peng-jiang@uiowa.edu