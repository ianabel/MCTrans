#!/usr/bin/env python3

from netCDF4 import Dataset
import numpy as np

import sys

file1 = ""
file2 = ""

if len(sys.argv) == 3:
    file1 = sys.argv[1]
    file2 = sys.argv[2]
else:
    print("Usage: ./nc_compare <netcdf filename> <netcdf filename>")
    sys.exit()

nc_root_1 = Dataset(file1, "r", format="NETCDF4")
nc_root_2 = Dataset(file2, "r", format="NETCDF4")

for variableName in nc_root_1.variables:
    v1 = nc_root_1.variables[variableName][...]
    v2 = nc_root_2.variables[variableName][...]
    if( not np.array_equal(v1,v2) ):
        sys.exit( 1 )

sys.exit( 0 )


