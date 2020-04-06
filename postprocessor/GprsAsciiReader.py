#!/usr/bin/env python
import pandas as pd
import numpy as np
import os

class GprsAsciiReader:
    """
    Reads gprs output in ascii format
    """
    # file_name = ""
    input_file = None
    n_blocks = 0                # for flow
    n_nodes = 0                 # for geomech
    current_time = 0.0
    data = pd.DataFrame()

    def __init__(self, file_name, n_mech_vert=0):
        """
        Need number of mech vertices as an arguments since ascii
        file does not specify that
        """
        self.input_file = open(file_name, "r")
        self.input_file.readline()         # ASCII VARIABLES header
        line = self.input_file.readline()  # number of blocks
        self.n_blocks = int(line.split()[-1]) # number of flow cells
        self.n_nodes = n_mech_vert            # number of vertices

    def __del__(self):
        if self.input_file is not None:
            self.input_file.close()

    def advanceTimeStep(self) -> bool:
        """
        returns true if was able to read a timestep
        returns false if eof.
        The data read is stored untill the next invocation
        """
        line = ""
        while (not line.strip()):             # skip empty
            line = self.input_file.readline() # Time  = ...
            if not line: return False

        self.current_time = float(line.split()[-1])
        line = self.input_file.readline() # table headers
        keys = line.split()
        # if we got geomechanics, it goes first
        if (keys[0] == "node"):
            assert self.n_nodes > 0, "Cannot yet inferm n geomech nodes"
            storage = np.zeros([self.n_nodes, len(keys)])
            for i in range(self.n_nodes):
                values = [float(x) for x in self.input_file.readline().split() ]
                storage[i, :] = values
            # skip till flow data
            line = self.input_file.readline() # Time  = ...
            while (not line.strip()):             # skip empty
                line = self.input_file.readline() # Time  = ...
                if not line: return False
            # read headers for flow
            line = self.input_file.readline()
            keys = line.split()

        if (keys[0] == 'cell'):   # read flow
            storage = np.zeros([self.n_blocks, len(keys)])
            for i in range(self.n_blocks):
                values = [float(x) for x in self.input_file.readline().split() ]
                storage[i, :] = values

        # put into dataframe
        self.data = pd.DataFrame(storage, columns=keys)
        return True

    def getTime(self) -> float:
        return self.current_time

    def getData(self) -> pd.DataFrame:
        """
        returns the data read by readTimeStep.
        """
        return self.data

    def getRelativePosition(self) -> float:
        current_pos = self.input_file.tell()
        size = os.stat(self.input_file.name)[6]
        return float(current_pos) / float(size)
