#!/bin/bash

for st in 10
do
for th in 100
do
  ./fsm4_s.exe /data/not_backed_up/shared/graphs/mico.lg 0.001 $st $th 1> output/fsm4_s_mico_001_"$st"_"$th".output 2> output/fsm4_s_mico_001_"$st"_"$th".err 
  ./fsm4_s.exe /data/not_backed_up/shared/graphs/mico.lg 0.005 $st $th 1> output/fsm4_s_mico_005_"$st"_"$th".output 2> output/fsm4_s_mico_005_"$st"_"$th".err
  ./fsm4_s.exe /data/not_backed_up/shared/graphs/mico.lg 0.01 $st $th 1> output/fsm4_s_mico_01_"$st"_"$th".output 2> output/fsm4_s_mico_01_"$st"_"$th".err
  ./fsm4_s.exe /data/not_backed_up/shared/graphs/mico.lg 0.05 $st $th 1> output/fsm4_s_mico_05_"$st"_"$th".output 2> output/fsm4_s_mico_05_"$st"_"$th".err
done
done
