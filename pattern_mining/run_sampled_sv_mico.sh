#!/bin/bash

for TH in 8
do
  echo $TH
  ./fsm4_s_sv.exe /data/not_backed_up/shared/graphs/mico.lg 0.001 $TH 1> output/fsm4_s_sv_mico_001_$TH.output 2> output/fsm4_s_sv_mico_001_$TH.err
  sleep 1

  ./fsm4_s_sv.exe /data/not_backed_up/shared/graphs/mico.lg 0.005 $TH 1> output/fsm4_s_sv_mico_005_$TH.output 2> output/fsm4_s_sv_mico_005_$TH.err
  sleep 1

  ./fsm4_s_sv.exe /data/not_backed_up/shared/graphs/mico.lg 0.01 $TH 1> output/fsm4_s_sv_mico_01_$TH.output 2> output/fsm4_s_sv_mico_01_$TH.err 
  sleep 1

  ./fsm4_s_sv.exe /data/not_backed_up/shared/graphs/mico.lg 0.05 $TH 1> output/fsm4_s_sv_mico_05_$TH.output 2> output/fsm4_s_sv_mico_05_$TH.err 
  sleep 1
done
