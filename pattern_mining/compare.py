g1 = {}
with open('tmp') as f:
    for l in f.readlines():
        t = l.split()
        t = [int(k) for k in t]
        p = sorted(t)
        if tuple(p) not in g1:
            g1[tuple(p)] = []
        g1[tuple(p)].append(t)

g2 = {}
with open('tmp1') as f:
    for l in f.readlines():
        t = l.split()
        t = [int(k) for k in t]
        p = sorted(t)
        if tuple(p) not in g2:
            g2[tuple(p)] = []
        g2[tuple(p)].append(t)

print (g1)
print (g2)

for x in g1:
    if x not in g2:
        print (x) 
