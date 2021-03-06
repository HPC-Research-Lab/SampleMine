#!/bin/bash

./fsm4_ns_np.exe /data/not_backed_up/shared/graphs/citeseer.lg 0.001 1> fsm_ns_np_citeseer_4_001.output 2> fsm_ns_np_citeseer_4_001.err 

./fsm4_ns_np.exe /data/not_backed_up/shared/graphs/citeseer.lg 0.005 1> fsm_ns_np_citeseer_4_005.output 2> fsm_ns_np_citeseer_4_005.err

./fsm4_ns_np.exe /data/not_backed_up/shared/graphs/citeseer.lg 0.01 1> fsm_ns_np_citeseer_4_01.output 2> fsm_ns_np_citeseer_4_01.err

./fsm4_ns_np.exe /data/not_backed_up/shared/graphs/citeseer.lg 0.05 1> fsm_ns_np_citeseer_4_05.output 2> fsm_ns_np_citeseer_4_05.err

echo "Done" | mail -s "Job1" peng-jiang@uiowa.edu

./fsm5_ns_np.exe /data/not_backed_up/shared/graphs/citeseer.lg 0.001 1> fsm_ns_np_citeseer_5_001.output 2> fsm_ns_np_citeseer_5_001.err 

./fsm5_ns_np.exe /data/not_backed_up/shared/graphs/citeseer.lg 0.005 1> fsm_ns_np_citeseer_5_005.output 2> fsm_ns_np_citeseer_5_005.err

./fsm5_ns_np.exe /data/not_backed_up/shared/graphs/citeseer.lg 0.01 1> fsm_ns_np_citeseer_5_01.output 2> fsm_ns_np_citeseer_5_01.err

./fsm5_ns_np.exe /data/not_backed_up/shared/graphs/citeseer.lg 0.05 1> fsm_ns_np_citeseer_5_05.output 2> fsm_ns_np_citeseer_5_05.err

echo "Done" | mail -s "Job2" peng-jiang@uiowa.edu

./fsm6_ns_np.exe /data/not_backed_up/shared/graphs/citeseer.lg 0.001 1> fsm_ns_np_citeseer_6_001.output 2> fsm_ns_np_citeseer_6_001.err 

./fsm6_ns_np.exe /data/not_backed_up/shared/graphs/citeseer.lg 0.005 1> fsm_ns_np_citeseer_6_005.output 2> fsm_ns_np_citeseer_6_005.err

./fsm6_ns_np.exe /data/not_backed_up/shared/graphs/citeseer.lg 0.01 1> fsm_ns_np_citeseer_6_01.output 2> fsm_ns_np_citeseer_6_01.err

./fsm6_ns_np.exe /data/not_backed_up/shared/graphs/citeseer.lg 0.05 1> fsm_ns_np_citeseer_6_05.output 2> fsm_ns_np_citeseer_6_05.err

echo "Done" | mail -s "Job3" peng-jiang@uiowa.edu
