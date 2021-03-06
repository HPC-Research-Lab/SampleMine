#!/bin/bash

for MT in 32  
do 
for ST in 4 16 36
do
  ./fsm5_s_ms.exe /data/not_backed_up/shared/graphs/com-orkut.lg $MT $ST 1> output/fsm5_s_ms_com-orkut_"$MT"_"$ST".output 2> output/fsm5_s_ms_com-orkut_"$MT"_"$ST".err
done
done
