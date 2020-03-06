#!/usr/bin/env python
import h5py
import pandas as pd
import numpy as np
import os

class GprsHDF5Reader:
    input_file = None
    n_blocks = 0
    current_step = 0
    data = pd.DataFrame()
    times = None

    def __init__(self, file_name):
        self.input_file = h5py.File(file_name, "r")
        self.times = self.input_file["RESTART"]["time"][:]
        self.n_blocks = self.input_file["FLOW_CELL"]["pres"].shape[0]

    def advanceTimeStep(self) -> bool:
        data = self.input_file["FLOW_CELL"]
        keys = data.keys()
        storage = np.zeros( [ self.n_blocks, len(keys) ] )
        for i, key in enumerate(keys):
            storage[:, i] = data[key][:, self.current_step]
        self.data = pd.DataFrame(storage, columns=keys)
        if self.current_step < len(self.times) - 1:
            self.current_step += 1
            return True
        else:
            return False

    def getTime(self) -> float:
        return self.times[ self.current_step ]

    def getData(self) -> pd.DataFrame:
        return self.data

    def getRelativePosition(self) -> float:
        return float(self.current_step) / float(len(self.times))

    def __del__(self):
        if self.input_file is not None:
            self.input_file.close()
