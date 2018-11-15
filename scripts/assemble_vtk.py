import os
import numpy as np

# case_path = "/home/ishovkun/sim/edfm-1frac/"
# case_path = "/home/ishovkun/sim/edfm-1frac-0/"
# case_path = "/home/ishovkun/sim/aquifer/"
# case_path = "/home/ishovkun/sim/edfm-dfm/"
case_path = "/home/ishovkun/sim/mech-tests/edfm-sneddon/"

res_mesh_file = case_path + "reservoir_mesh.vtk"
edfm_mesh_file = case_path + "efrac.vtk"
dfm_mesh_file = case_path + "dfm.vtk"
results_file = case_path + "OUTPUT.vars.txt"
vtk_dir = case_path + "OUTPUT.vtk_output/"
gm_sda_file = case_path + "gm_SDA.txt"
output_dir = case_path + "field"

n_dfm = 0
try:
    with open(dfm_mesh_file, "r") as f:
        line = f.readline().rstrip()
        split = line.split()
        while (split[0] != "CELL_TYPES"):
            line = f.readline().rstrip()
            split = line.split()
            if (len(split) == 0):
                split.append(-1)
                continue
        n_dfm = int(split[1])
except FileNotFoundError:
    pass

print("n dfm segments = ", n_dfm)

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
                # print(data)
                n_elements = int(data[1])
                for i in range(n_elements):
                    line = f.readline().rstrip()
                    value = int(line.split()[0])
                    if (value == 12 or value == 13):
                        n_cells += 1
                    else:
                        print("unknows value ", value)
                        exit(0)
        counter += 1

print("n_cells =", n_cells)

n_edfm = 0
try:
    with open(edfm_mesh_file, "r") as f:
        line = f.readline().rstrip()
        split = line.split()
        while (split[0] != "CELL_TYPES"):
            line = f.readline().rstrip()
            split = line.split()
            if (len(split) == 0):
                split.append(-1)
                continue
        n_edfm = int(split[1])
except FileNotFoundError:
    pass

print("n edfm = ", n_edfm)

print("reading ASCII flow data")
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

        # print(split)
        if (split[0] == "node"):  # skip mechanics
            while(len(split) != 0):  # data loop within time step
                line = f.readline().rstrip()
                split = line.split()
                # print(split)
            line = f.readline().rstrip()
            times.pop()
            continue

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

        all_data.append(data)
        line = f.readline().rstrip()

# get edfm jumps from vtk files and add them to data
# sda_data = np.zeros([n_edfm, 2])
sda_data  = []
if (n_edfm > 0):
    with open (gm_sda_file, "r") as f:
        found = False
        n_sda_cells = 0
        counter = 0
        sda_cells = []
        for word in [word for line in f for word in line.split()]:
            if (word == "GM_EFRAC_CELLS"):
                found = True
                continue
            if (word == "/"):
                break
            if (n_sda_cells == 0):
                n_sda_cells = int(word)
                continue
            if (n_cells > 0):
                sda_cells.append(int(word) - 1)
                counter += 1
                if (counter == n_sda_cells):
                    counter = 0


    print(sda_cells)

    file_counter = 0
    # print(sorted(os.listdir(vtk_dir)))
    for fname in sorted(os.listdir(vtk_dir)):
        sda_data.append(np.zeros([n_edfm, 2]))
        if ("block_scalars" in fname):  # found vtk block file
            print ("Reading ", fname)
            with open(vtk_dir + fname, "r") as f:
                # edfm_cell_index = n_cells
                cell_counter = 0
                sda_counter = 0
                found_segment = False
                skip_line = True
                for line in f:
                    split = line.split()
                    if (len(split) == 0): continue

                    if (len(split) > 1):
                        if(split[1] == "Jump_n"):
                            found_segment = True
                            continue

                    if (found_segment):
                        if (skip_line):  # skip LOOKUP_TABLE HSV
                            skip_line = False
                            continue
                        else:   # start counting cells
                            # print(split)
                            if (cell_counter in sda_cells):
                                sda_data[file_counter][sda_counter, 0] = float(split[0])
                                sda_counter += 1

                            cell_counter += 1
                            if (cell_counter == n_cells):
                                found_segment = False
                                cell_counter = 0
                                skip_line = True

print("n_volumes = ",                n_volumes,
      "|| n_dfm + n_edfm + n_cells = ", n_dfm + n_edfm + n_cells)
if (n_volumes != n_dfm + n_edfm + n_cells):
    print("bug # volumes")
    exit(0)

print("reading reservoir geometry")
res_msh_txt = None
with open(res_mesh_file, "r") as f:
    res_msh_txt = f.read()

edfm_msh_txt = None
if (n_edfm > 0):
    print("reading edfm geometry")
    with open(edfm_mesh_file, "r") as f:
        edfm_msh_txt = f.read()

dfm_msh_txt = None
if (n_dfm > 0):
    print("reading dfm geometry")
    with open(dfm_mesh_file, "r") as f:
        dfm_msh_txt = f.read()


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



print("saving reservoir data")
# save reservoir vtk
for i in range(len(times)):
    with open(output_dir + "/" + "blocks-" + str(i) + ".vtk", "w") as f:
        f.write(res_msh_txt)
        f.write("CELL_DATA " + str(n_cells) + "\n")
        for j in range(len(headers[i])):
            header = headers[i][j].replace("=", "_")
            f.write("SCALARS\t" + header + "\tfloat\n")
            f.write("LOOKUP_TABLE HSV\n")
            # print(all_data[i].shape)
            # for k in range(n_frac_elements):
            for k in range(n_cells):
                f.write("%.5f\n" % all_data[i][n_dfm + k, j])

if (n_dfm > 0):
    print("saving DFM data")
    for i in range(len(times)):
        with open(output_dir + "/" + "dfm-" + str(i) + ".vtk", "w") as f:
            f.write(dfm_msh_txt)
            f.write("CELL_DATA " + str(n_dfm) + "\n")
            for j in range(len(headers[i])):
                header = headers[i][j].replace("=", "_")
                f.write("SCALARS\t" + header + "\tfloat\n")
                f.write("LOOKUP_TABLE HSV\n")
                for k in range(n_dfm):
                    f.write("%.5f\n" % all_data[i][k, j])

if (n_edfm > 0):
    print("saving EDFM data")
    for i in range(len(times)):
        with open(output_dir + "/" + "edfm-" + str(i) + ".vtk", "w") as f:
            f.write(edfm_msh_txt)
            f.write("CELL_DATA " + str(n_edfm) + "\n")
            for j in range(len(headers[i])):
                header = headers[i][j].replace("=", "_")
                f.write("SCALARS\t" + header + "\tfloat\n")
                f.write("LOOKUP_TABLE HSV\n")
                for k in range(n_edfm):
                    f.write("%.5f\n" % all_data[i][n_dfm + n_cells + k, j])

            # write sda data
            for j in range(1):
                header = "jump_n"
                f.write("SCALARS\t" + header + "\tfloat\n")
                f.write("LOOKUP_TABLE HSV\n")
                for k in range(n_edfm):
                    f.write("%.5f\n" % sda_data[i][k, j])
