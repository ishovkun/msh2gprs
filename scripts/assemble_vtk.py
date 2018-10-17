import os
import numpy as np

# case_path = "/home/production/sim/edfm-2frac/"
case_path = "/home/production/sim/aquifer/"
res_mesh_file = case_path + "reservoir_mesh.vtk"
frac_mesh_file = case_path + "efrac.vtk"
results_file = case_path + "OUTPUT.vars.txt"
output_dir = case_path + "field"
# res_mesh_file = "/home/ishovkun/sim/edfm-2frac/reservoir_mesh.vtk"
# frac_mesh_file = "/home/ishovkun/sim/edfm-2frac/efrac.vtk"
# results_file = "/home/ishovkun/sim/edfm-2frac/OUTPUT.vars.txt"
# output_dir = "/home/ishovkun/sim/edfm-2frac/field"

n_cells = 0
with open(res_mesh_file, "r") as f:
    break_flag = False
    counter = 0
    line = True
    # while (not break_flag or counter > 30):
    # while (line and counter < 30):
    while (line):               # outer loop
        line = f.readline()
        data = line.split()
        if (len(data) > 0):
            if (data[0] == "CELL_TYPES"):
                print(data)
                # n_cells = int(data[1])
                n_elements = int(data[1])
                for i in range(n_elements):
                    line = f.readline().rstrip()
                    value = int(line.split()[0])
                    if (value == 12):
                        n_cells += 1
        counter += 1

print("n_cells =", n_cells)
# exit(0)
# print("quit")

times = []
headers = []
all_data = []
n_volumes = 0
with open(results_file, "r") as f:
    f.readline()                # ASCII VARIABLES header
    line = f.readline()             # number of blocks
    n_volumes = int(line.split()[-1])
    line = f.readline().rstrip()
    while (line):               # time loop
        split = line.split()
        t = float(split[-1])
        times.append(t)

        line = f.readline().rstrip()
        split = line.split()
        n_cols = len(split) - 1
        data = np.zeros([n_volumes, n_cols])
        headers.append(split[1:])

        # print(split)
        row = 0
        while(len(split) != 0):  # data loop within time step
            line = f.readline().rstrip()
            split = line.split()
            if (len(split) == 0):
                continue

            for i in range(1, n_cols+1):
                data[row, i-1] = float(split[i])
                # print(split[i])
            row += 1

        # print(data)
        # exit(0)
        all_data.append(data)
        line = f.readline().rstrip()

n_frac_elements = n_volumes - n_cells

frac_msh_txt = None
with open(frac_mesh_file, "r") as f:
    frac_msh_txt = f.read()

res_msh_txt = None
with open(res_mesh_file, "r") as f:
    res_msh_txt = f.read()


# save efrac vtk
import shutil
try:
    shutil.rmtree(output_dir)
except:
    pass

try:
    os.mkdir(output_dir)
except FileExistsError:
    pass

for i in range(len(times)):
    with open(output_dir + "/" + "frac-" + str(i) + ".vtk", "w") as f:
        f.write(frac_msh_txt)
        f.write("CELL_DATA " + str(n_frac_elements) + "\n")
        for j in range(len(headers[i])):
            header = headers[i][j].replace("=", "_")
            f.write("SCALARS\t" + header + "\tfloat\n")
            f.write("LOOKUP_TABLE HSV\n")
            print(all_data[i].shape)
            for k in range(n_frac_elements):
                f.write("%.5f\n" % all_data[i][k+n_cells, j])


# save reservoir vtk
for i in range(len(times)):
    with open(output_dir + "/" + "blocks-" + str(i) + ".vtk", "w") as f:
        f.write(res_msh_txt)
        # f.write("CELL_DATA " + str(n_frac_elements) + "\n")
        f.write("CELL_DATA " + str(n_cells) + "\n")
        for j in range(len(headers[i])):
            header = headers[i][j].replace("=", "_")
            f.write("SCALARS\t" + header + "\tfloat\n")
            f.write("LOOKUP_TABLE HSV\n")
            print(all_data[i].shape)
            # for k in range(n_frac_elements):
            for k in range(n_cells):
                f.write("%.5f\n" % all_data[i][k, j])
