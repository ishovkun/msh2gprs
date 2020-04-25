from hybmeshpack import hmscript as hm
import numpy as np
# This script was created for hybmesh 0.4.5.
# If running version of hybmesh is not compatible with this version
# an exception will be raised at this line.
hm.check_compatibility("0.4.5", 2)

# ====================== Input data
# define wells coordinates

# area and fracture coordinates
area_pts = [-0.5, -0.5], [0.5, 0.5]

# meshing options
npoint = 4
dom_size = 1. / npoint   # grid cell size in outer domain

area = hm.add_rect_contour(*area_pts)


# ======================= FEM mesh
parea = hm.partition_contour(area, 'const', dom_size)
fin5 = hm.pebi_fill(parea)
hm.export_grid_vtk(fin5, "grid.vtk")
zcoords = [0, 2.5, 5]
#  zcoords = np.linspace(0, 5, 5*npoint)
fin7 = hm.extrude_grid(fin5, zcoords)

hm.export3d_grid_vtk(fin7, "grid.vtk")
