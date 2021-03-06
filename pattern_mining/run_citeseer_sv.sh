#!/bin/bash


./fsm7_ns_sv.exe /data/not_backed_up/shared/graphs/citeseer.lg 0.001 1> fsm_ns_sv_citeseer_7_001.output 2> fsm_ns_sv_citeseer_7_001.err 

./fsm7_ns_sv.exe /data/not_backed_up/shared/graphs/citeseer.lg 0.005 1> fsm_ns_sv_citeseer_7_005.output 2> fsm_ns_sv_citeseer_7_005.err

./fsm7_ns_sv.exe /data/not_backed_up/shared/graphs/citeseer.lg 0.01 1> fsm_ns_sv_citeseer_7_01.output 2> fsm_ns_sv_citeseer_7_01.err

./fsm7_ns_sv.exe /data/not_backed_up/shared/graphs/citeseer.lg 0.05 1> fsm_ns_sv_citeseer_7_05.output 2> fsm_ns_sv_citeseer_7_05.err

echo "Done" | mail -s "Job" peng-jiang@uiowa.edu
