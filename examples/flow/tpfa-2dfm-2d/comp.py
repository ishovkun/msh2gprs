#!/usr/bin/env python
import numpy as np

fname1 = "fl_face_data.txt"
fname2 = "fl_face_data-old.txt"
# fname1 = "../tpfa-nofrac-2d/fl_face_data_old.txt"
# fname2 = "../tpfa-nofrac-2d/fl_face_data.txt"

def read_con_data(fname):
    data = []
    with open (fname, "r") as f:
        n_conns = 0
        for i, line in enumerate(f):
            if (i == 0): continue
            elif (i == 1):
                n_conns = int(line)
                # data = np.zeros([n_conns, 3], names=["i", "j", "T"])
            elif (i == n_conns + 2): break
            else:
                nums = [float(x) for x in line.split()]
                # data[i - 2, :] = nums
                data.append([int(nums[0]), int(nums[1]), nums[2]])
    return np.array(data)

def print_conn(d):
    for i in range(d.shape[0]):
        print ("%d\t\t%d\t\t%f" % ((int(d[i,0]), int(d[i,1]), d[i,2]  )))

def sort_by_second(d):
    start = 0
    end = 0
    n_rows = d.shape[0]
    while (start < n_rows - 1):
        # deternime sorting bounds
        for i in range(start+1, n_rows):
            end = i
            if (d[i, 0] == d[i-1, 0]): continue
            else: break
        # sort
        order = d[start:end,1].argsort(kind='stable')
        d[start:end, :] = d[start:end, :][order]
        start = end
    return d


d1 = read_con_data(fname1)
d1 = d1[d1[:,0].argsort(kind='stable')]
sort_by_second(d1)
print("old"); print_conn(d1)

d2 = read_con_data(fname2)
d2 = d2[d2[:,0].argsort(kind='stable')]
sort_by_second(d2)
print("new"); print_conn(d2)

for i in range(d1.shape[0]):
    print( "%d %d %1.3f "%(d1[i, 0] , d1[i, 1],
                           abs(d1[i, -1] - d2[i, -1]) / d2[i, -1] * 100.))
