#!/bin/bash

./fsm_ns.exe /data/not_backed_up/shared/graphs/citeseer.lg 0.001 1> fsm_ns_citeseer_6_001.output 2> fsm_ns_citeseer_6_001.err 

./fsm_ns.exe /data/not_backed_up/shared/graphs/citeseer.lg 0.005 1> fsm_ns_citeseer_6_005.output 2> fsm_ns_citeseer_6_005.err

./fsm_ns.exe /data/not_backed_up/shared/graphs/citeseer.lg 0.01 1> fsm_ns_citeseer_6_01.output 2> fsm_ns_citeseer_6_01.err

./fsm_ns.exe /data/not_backed_up/shared/graphs/citeseer.lg 0.05 1> fsm_ns_citeseer_6_05.output 2> fsm_ns_citeseer_6_05.err
