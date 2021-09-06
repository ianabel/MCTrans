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

ti_var = nc_root.variables["IonTemperature"]
te_var = nc_root.variables["ElectronTemperature"]
t_var = nc_root.variables["t"]

plt.figure(figsize=(10,6),dpi=210)
plt.title("Temperatures" )
plt.xlabel("t (seconds)")
plt.ylabel("Temperature (keV)" )
plt.plot(t_var[:],ti_var[:],label='T_i',linewidth=1)
plt.plot(t_var[:],te_var[:],label='T_e',linewidth=1)
plt.legend()
plt.show()


j_var = nc_root.variables["Current"]
v_var = nc_root.variables["Voltage"]

fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('time (s)')
ax1.set_ylabel('Voltage (V)', color=color)
ax1.plot(t_var, v_var, label='V', color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('Radial Current (A)', color=color)  # we already handled the x-label with ax1
ax2.plot(t_var, j_var, color=color)
ax2.tick_params(axis='y', labelcolor=color)

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.legend()
plt.show()

