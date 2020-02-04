#!/usr/bin/env python
from GprsAsciiReader import GprsAsciiReader
import yaml
import vtk
import os, shutil
from vtk.util import numpy_support
from printProgressBar import printProgressBar
# import pyvtk

class Postprocessor:
    """
    Main class of the program
    """
    def __init__(self, case_path):
        self.case_path = case_path
        self.config_file_name = "postprocessor_config.yaml"
        self.config = self.readConfig_()
        self.matrix_flow_grid_reader = vtk.vtkDataSetReader()
        self.edfm_flow_grid_reader = vtk.vtkDataSetReader()
        self.dfm_flow_grid_reader = vtk.vtkDataSetReader()
        self.writer = vtk.vtkDataSetWriter()
        self.readGeometry_()
        self.output_file_number = 0
        self.output_dir = self.case_path + self.config["output_directory"]

    def run(self):
        gprs_reader = GprsAsciiReader(self.case_path + "OUTPUT.vars.txt")
        self.prepareOutputDirectory_()

        printProgressBar(0, 1, prefix = 'Progress:', suffix = 'Complete', length = 20)
        while gprs_reader.readTimeStep():
           t = gprs_reader.getTime()
           data = gprs_reader.getData()
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
        # dfm grid
        if len(self.config["dfm_cell_to_flow_dof"]) > 0:
            reader = self.dfm_flow_grid_reader
            vtk_file_path = self.case_path + self.config["dfm_flow_grid_file"]
            reader.SetFileName(vtk_file_path)
            reader.Update()
            assert reader.GetOutput().GetNumberOfCells() == len(self.config["dfm_cell_to_flow_dof"])
        # edfm grid
        if len(self.config["edfm_cell_to_flow_dof"]) > 0:
            reader = self.edfm_flow_grid_reader
            vtk_file_path = self.case_path + self.config["edfm_grid_file"]
            reader.SetFileName(vtk_file_path)
            reader.Update()
            assert reader.GetOutput().GetNumberOfCells() == len(self.config["edfm_cell_to_flow_dof"])

    def saveReservoirData_(self, t, data):
        assert data.shape[0] == (len(self.config["matrix_cell_to_flow_dof"]) +
                                 len(self.config["dfm_cell_to_flow_dof"]) +
                                 len(self.config["edfm_cell_to_flow_dof"])), "Data size " + \
                                 str(data.shape[0]) + " != " +\
                        str(len( self.config["matrix_cell_to_flow_dof"] ) )+ " " + \
                        str(len( self.config["dfm_cell_to_flow_dof"] )) + " " +\
                        str(len(self.config["edfm_cell_to_flow_dof"]))

        # Extract flor, edfm, and dfm data
        self.addDataToReader_(self.matrix_flow_grid_reader, data, self.config["matrix_cell_to_flow_dof"])
        self.addDataToReader_(self.edfm_flow_grid_reader, data, self.config["edfm_cell_to_flow_dof"])
        self.addDataToReader_(self.dfm_flow_grid_reader, data, self.config["dfm_cell_to_flow_dof"])
        # save output
        self.writeFile_(self.matrix_flow_grid_reader, "matrix-%d"%self.output_file_number)
        self.writeFile_(self.edfm_flow_grid_reader, "edfm-%d"%self.output_file_number)
        self.writeFile_(self.dfm_flow_grid_reader, "dfm-%d"%self.output_file_number)
        self.output_file_number += 1

    def addDataToReader_(self, reader, data, mapping):
        output = reader.GetOutput()
        for key in data.keys():
            x = numpy_support.numpy_to_vtk(data[key].values[mapping])
            x.SetName(key)
            output.GetCellData().AddArray(x)

    def writeFile_(self, reader, file_name):
        self.writer.SetFileName(self.output_dir + "/" + file_name + ".vtk")
        self.writer.SetInputData(reader.GetOutput())
        self.writer.Write()

    def prepareOutputDirectory_(self):
        if (os.path.isdir(self.output_dir)):
            shutil.rmtree(self.output_dir)
        os.mkdir(self.output_dir)

if __name__ == "__main__":
    case_path = "/home/ishovkun/sim/embedded_fractures/hybrid/single_dfm_single_edfm/"
    proc = Postprocessor(case_path)
    proc.run()
