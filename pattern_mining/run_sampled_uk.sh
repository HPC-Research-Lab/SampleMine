#!/bin/bash

for TH in 2 
do
  echo $TH
  ./fsm5_s_ms.exe /data/not_backed_up/shared/graphs/uk-2005/uk-2005.lg $TH 1> output/fsm5_s_ms_uk-2005_$TH.output 2> output/fsm5_s_ms_uk-2005_$TH.err
done
