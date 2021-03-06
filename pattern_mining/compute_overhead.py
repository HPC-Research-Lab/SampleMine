import sys
import re

r_overhead = 0
a_overhead = 0

with open(sys.argv[1]) as f:
  for l in f.readlines():
    m = re.findall(r'merging', l)
    if m:
      num = re.findall(r'[0-9]+', l)
      r_overhead += int(num[0]) * int(num[1])
      a_overhead += min(int(num[0]), int(sys.argv[2])) * min(int(num[1]), int(sys.argv[2])) 

print (r_overhead / a_overhead)
