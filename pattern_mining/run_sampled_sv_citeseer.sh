#!/bin/bash

for TH in 2 
do
  echo $TH
  #./fsm9_s_sv.exe /data/not_backed_up/shared/graphs/citeseer.lg 0.001 $TH 1> output/fsm9_s_sv_citeseer_001_$TH.output 2> output/fsm9_s_sv_citeseer_001_$TH.err
#  sleep 1

  ./fsm9_s_sv.exe /data/not_backed_up/shared/graphs/citeseer.lg 0.005 $TH 1> output/fsm9_s_sv_citeseer_005_$TH.output 2> output/fsm9_s_sv_citeseer_005_$TH.err
  sleep 1

 ./fsm9_s_sv.exe /data/not_backed_up/shared/graphs/citeseer.lg 0.01 $TH 1> output/fsm9_s_sv_citeseer_01_$TH.output 2> output/fsm9_s_sv_citeseer_01_$TH.err 
   sleep 1

  ./fsm9_s_sv.exe /data/not_backed_up/shared/graphs/citeseer.lg 0.05 $TH 1> output/fsm9_s_sv_citeseer_05_$TH.output 2> output/fsm9_s_sv_citeseer_05_$TH.err 
  sleep 1
done
