#!/usr/bin/env python
from GprsAsciiReader import GprsAsciiReader
from GprsHDF5Reader import GprsHDF5Reader
import yaml
import vtk
import os, shutil
from vtk.util import numpy_support
from printProgressBar import printProgressBar
import numpy as np
import pandas as pd

# import pyvtk

class Postprocessor:
    """
    Main class of the program
    """
    def __init__(self, case_path):
        self.case_path = case_path
        self.config_file_name = "postprocessor_config.yaml"
        self.matrix_flow_grid_reader = vtk.vtkDataSetReader()
        self.dfm_flow_grid_reader = vtk.vtkDataSetReader()
        self.matrix_mech_grid_reader = vtk.vtkDataSetReader()
        self.matrix_mech_vtk_ouput_reader = vtk.vtkDataSetReader()
        self.writer = vtk.vtkDataSetWriter()
        self.output_file_number = 0
        self.config = self.readConfig_()
        self.output_dir = self.case_path + self.config["output_directory"]
        self.readGeometry_()
        self.has_mechanics = self.checkMechanicsVTKOutput()

    def run(self):
        gprs_reader = self.makeReader_()
        self.prepareOutputDirectory_()

        # printProgressBar(0, 1, prefix = 'Progress:', suffix = 'Complete', length = 20)
        while gprs_reader.advanceTimeStep():
           t = gprs_reader.getTime()
           data = gprs_reader.getData()

           if (self.has_mechanics):
               mech_data = gprs_reader.getMechData()
               vtk_mech_matrix_vtk_data = self.readMechVTKData("block_scalars")
               vtk_mech_fracture_data = None
               if (self.config["dfm_cell_to_flow_dof"]):
                   vtk_mech_fracture_data = self.readMechVTKData("frac_scalars")
               self.saveReservoirData_(t, data, mech_data=mech_data,
                                       mech_vtk_data=vtk_mech_matrix_vtk_data,
                                       mech_vtk_frac_data=vtk_mech_fracture_data)
           else:
               self.saveReservoirData_(t, data)

           printProgressBar(gprs_reader.getRelativePosition(), 1, prefix = 'Progress:', suffix = 'Complete', length = 20)
        printProgressBar(1, 1, prefix = 'Progress:', suffix = 'Complete', length = 20)

    def readConfig_(self):
        if (not os.path.isdir(self.case_path)):
            raise IsADirectoryError("%s does not exist" % self.case_path)
        with open(self.case_path + self.config_file_name, "r") as f:
            return yaml.load(f, Loader=yaml.CLoader)
        if (self.config is None):
           raise FileExistsError(self.case_path + self.config_file_name)

    def readGeometry_(self):
        # reservoir grid
        reader = self.matrix_flow_grid_reader
        vtk_file_path = self.case_path + self.config["flow_reservoir_grid_file"]
        reader.SetFileName(vtk_file_path)
        reader.Update()
        assert reader.GetOutput().GetNumberOfCells() == len(self.config["matrix_cell_to_flow_dof"])
        n_cells_flow_grid = reader.GetOutput().GetNumberOfCells()
        # dfm grid
        if len(self.config["dfm_cell_to_flow_dof"]) > 0:
            reader = self.dfm_flow_grid_reader
            vtk_file_path = self.case_path + self.config["fracture_grid_file"]
            reader.SetFileName(vtk_file_path)
            reader.Update()
            assert reader.GetOutput().GetNumberOfCells() == len(self.config["dfm_cell_to_flow_dof"])
        # geomechanics
        reader = self.matrix_mech_grid_reader
        vtk_file_path = self.case_path + self.config['mech_reservoir_grid_file']
        reader.SetFileName(vtk_file_path)
        reader.Update()
        self.separate_mech_file = False
        self.n_mech_vertices = reader.GetOutput().GetNumberOfPoints()
        n_cells_mech_grid = reader.GetOutput().GetNumberOfCells()
        if (n_cells_flow_grid != n_cells_mech_grid):
            self.separate_mech_file = True

    def saveReservoirData_(self, t, data, mech_data=None, mech_vtk_data=None,
                           mech_vtk_frac_data=None):
        assert data.shape[0] == (len(self.config["matrix_cell_to_flow_dof"]) +
                                 len(self.config["dfm_cell_to_flow_dof"])), \
                                 "Data size " + \
                                 str(data.shape[0]) + " != " +\
                        str(len( self.config["matrix_cell_to_flow_dof"] ) )+ " " + \
                        str(len( self.config["dfm_cell_to_flow_dof"] )) + " " +\
                        str(len(self.config["edfm_cell_to_flow_dof"]))

        # Extract flor, edfm, and dfm data
        if (not self.separate_mech_file):
            flow_reader = self.matrix_mech_grid_reader
        else:
            flow_reader = self.matrix_flow_grid_reader

        self.addDataToReader_(flow_reader, data, self.config["matrix_cell_to_flow_dof"], t)
        # flow_reader
        if len(self.config["dfm_cell_to_flow_dof"]) > 0:
            self.addDataToReader_(self.dfm_flow_grid_reader, data,
                                  self.config["dfm_cell_to_flow_dof"], t)
        if (mech_data is not None):
            self.addMechDataToReader(self.matrix_mech_grid_reader, mech_vtk_data, "cell")
            self.addMechDataToReader(self.matrix_mech_grid_reader, mech_data, "point")
            if ( mech_vtk_frac_data is not None ):
                self.addDataToReader_(self.dfm_flow_grid_reader, mech_vtk_frac_data,
                                      np.r_[0:len(mech_vtk_frac_data)], t)

        # save output
        if (not self.separate_mech_file):
            self.writeFile_(self.matrix_mech_grid_reader, "matrix-%d"%self.output_file_number)
        else:
            self.writeFile_(self.matrix_mech_grid_reader, "matrix-mech-%d"%self.output_file_number)
            self.writeFile_(self.matrix_flow_grid_reader, "matrix-%d"%self.output_file_number)

        if len(self.config["edfm_cell_to_flow_dof"]) > 0:
            self.writeFile_(self.edfm_flow_grid_reader, "edfm-%d"%self.output_file_number)
        if len(self.config["dfm_cell_to_flow_dof"]) > 0:
            self.writeFile_(self.dfm_flow_grid_reader, "dfm-%d"%self.output_file_number)
        self.output_file_number += 1

    def addDataToReader_(self, reader, data, mapping, time):
        output = reader.GetOutput()
        for key in data.keys():
            x = numpy_support.numpy_to_vtk(data[key].values[mapping])
            x.SetName(key)
            output.GetCellData().AddArray(x)
        # make time array the size of the n_cells; otherwise,
        # I cannot plot selection vs real time
        time_arr = numpy_support.numpy_to_vtk(time*np.ones(len( mapping )))
        time_arr.SetName("Time_")
        output.GetCellData().AddArray(time_arr)
        # output.GetFieldData().AddArray(time_arr)


    def addMechDataToReader(self, reader, data, data_type="point"):
        output = reader.GetOutput()
        for key in data.keys():
            x = numpy_support.numpy_to_vtk(data[key].values)
            x.SetName(key)
            if (data_type=="point"):
                output.GetPointData().AddArray(x)
            else:
                output.GetCellData().AddArray(x)

    def writeFile_(self, reader, file_name):
        self.writer.SetFileName(self.output_dir + "/" + file_name + ".vtk")
        self.writer.SetInputData(reader.GetOutput())
        self.writer.Write()

    def prepareOutputDirectory_(self):
        if (os.path.isdir(self.output_dir)):
            shutil.rmtree(self.output_dir)
        os.mkdir(self.output_dir)

    def makeReader_(self):
        if (os.path.isfile(self.case_path + "OUTPUT.vars.h5")):
            return GprsHDF5Reader(self.case_path + "OUTPUT.vars.h5")
        elif (os.path.isfile(self.case_path + "OUTPUT.vars.txt")):
            return GprsAsciiReader(self.case_path + "OUTPUT.vars.txt", self.n_mech_vertices)
        else:
           raise FileExistsError("Could not find GPRS output file")

    def checkMechanicsVTKOutput(self):
        '''
        Checks if mech vtk output directory exists
        '''
        # if (not os.path.isdir(self.case_path)):
        if os.path.isdir( self.case_path + "/OUTPUT.vtk_output" ):
            print("Mechanics vtk file found")
            return True
        else:
            print("No mechanics since vtk file found")
            return False

    def readMechVTKData(self, fname):
        reader = self.matrix_mech_vtk_ouput_reader
        if (self.output_file_number == 0):
            fnum_str = "%09d"%(self.output_file_number)
            # fnum_str = "%09d"%(self.output_file_number + 1)
        else:
            fnum_str = "%09d"%(self.output_file_number - 1)

        vtk_file_path = self.case_path + "/OUTPUT.vtk_output/" + fname + \
                        "." + fnum_str + ".vtk"
        # print("reading ", vtk_file_path)
        reader.SetFileName(vtk_file_path)
        reader.ReadAllScalarsOn()       # without this vtk reads only one array
        reader.Update()
        data = reader.GetOutput().GetCellData()
        n_vtk_arrays = data.GetNumberOfArrays()
        n_cells = reader.GetOutput().GetNumberOfCells()
        storage = np.zeros([n_cells, n_vtk_arrays])
        keys = []
        for i in range(n_vtk_arrays):
            keys.append(data.GetArrayName(i))
            storage[:, i] = numpy_support.vtk_to_numpy(data.GetArray(i))

        return pd.DataFrame(storage, columns=keys)


if __name__ == "__main__":
    case_path = "/home/ishovkun/sim/embedded_fractures/hybrid/3frac-dfm/"
    proc = Postprocessor(case_path)
    proc.run()
