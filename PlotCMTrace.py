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
    print("Usage: ./PlotMCTrace.py <netcdf filename>")
    sys.exit()

nc_root = Dataset(filename, "r", format="NETCDF4")

voltages_var = nc_root.variables["Voltage"]
t_var = nc_root.variables["t"]

print(t_var[10:].size)

fig, (ax1,ax3) = plt.subplots(2,1,sharex=True,figsize=(10,6),dpi=210)

fig.suptitle("Time Traces" )
ax1.set_xlabel("t / seconds")
ax1.set_ylabel("Volts" )
ax1.tick_params(axis='y', labelcolor='tab:blue')

ax1.plot(t_var[10:],voltages_var[10:],label = 'Applied Voltage',color = 'tab:red')

ax1.legend(loc='upper left')
ax2 = ax1.twinx()

ti_var = nc_root.variables["IonTemperature"];
te_var = nc_root.variables["ElectronTemperature"];


ax2.set_ylabel("keV")
ax2.tick_params(axis='y', labelcolor='tab:blue')
ax2.plot(t_var[10:],ti_var[10:],label = 'Ion Temperature',color = 'tab:blue')
ax2.plot(t_var[10:],te_var[10:],label = 'Electron Temperature',color = 'tab:green')
ax2.legend(loc='lower right')

fig.tight_layout()

ax3.set_xlabel("t / seconds")
ax3.set_ylabel("Volts" )
ax3.tick_params(axis='y', labelcolor='tab:blue')

ax3.plot(t_var[10:],voltages_var[10:],label = 'Applied Voltage',color = 'tab:red')

j_var = nc_root.variables["Current"]
ax4 = ax3.twinx()
ax4.set_ylabel("A")
ax4.tick_params(axis='y', labelcolor='tab:blue')
ax4.plot(t_var[10:],j_var[10:],label = 'Current',color = 'tab:blue')
ax4.legend(loc='lower right')




plt.show()

