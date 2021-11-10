#!/usr/bin/python3

from netCDF4 import Dataset
import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import sys

filename = ""

if len(sys.argv) == 2:
    filename = sys.argv[1]
else:
    print("Usage: ./plot_time_traces.py <netcdf filename>")
    sys.exit()

nc_root = Dataset(filename, "r", format="NETCDF4")

volt_var = nc_root.variables["Voltage"]
t_var = nc_root.variables["Time"]

plt.figure(figsize=(10,6),dpi=210)
plt.title("Imposed Voltage" )
plt.xlabel("t (seconds)")
plt.ylabel("Voltage (V)" )
plt.plot(t_var[:],volt_var[:],label='V',linewidth=1)
plt.legend()
plt.show()
