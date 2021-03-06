#!/bin/bash

for TH in 1 2 4
do
  echo $TH
  ./fsm5_s_sv.exe /data/not_backed_up/shared/graphs/com-orkut.lg $TH 1> output/fsm5_s_sv_com-orkut_$TH.output 2> output/fsm5_s_sv_com-orkut_$TH.err
done
