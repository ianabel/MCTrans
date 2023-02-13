#!/usr/bin/env python3

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

savefile = filename.split('/')[-1].split('.')[0]

nc_root = Dataset(filename, "r", format="NETCDF4")

voltages_var = nc_root.variables["Voltage"]
t_var = nc_root.variables["t"]
t_var = t_var[:]*1e3

t_start_idx = 5
figsize = (8,8)

print(t_var[t_start_idx:].size)

fig, axs = plt.subplots(3, 1, sharex=True, figsize=figsize,dpi=250)

axs[0].set_ylabel("kV")
axs[0].tick_params(axis='y', labelcolor='tab:red')

axs[0].plot(t_var[t_start_idx:],voltages_var[t_start_idx:]/1000,label = 'Applied Voltage',color = 'tab:red')

axs[0].legend(loc='lower left')
ax2 = axs[0].twinx()

ti_var = nc_root.variables["IonTemperature"]
te_var = nc_root.variables["ElectronTemperature"]

ax2.set_ylabel("keV")
ax2.tick_params(axis='y', labelcolor='tab:blue')
ax2.plot(t_var[t_start_idx:],ti_var[t_start_idx:],label = 'Ion Temperature',color = 'tab:blue')
ax2.plot(t_var[t_start_idx:],te_var[t_start_idx:],label = 'Electron Temperature',color = 'tab:green')
ax2.legend(loc='upper right')

axs[1].set_ylabel("W/m$^3$" )
axs[1].tick_params(axis='y', labelcolor='tab:blue')

viscous_var = nc_root.variables["ViscousTorque"]
axs[1].plot(t_var[t_start_idx:],viscous_var[t_start_idx:],label = 'Viscous Torque')
parallel_var = nc_root.variables["ParAngMomLoss"]
axs[1].plot(t_var[t_start_idx:],parallel_var[t_start_idx:],label = 'Parallel Angular Momentum Loss')
cx_var = nc_root.variables["CXAngMomLoss"]
axs[1].plot(t_var[t_start_idx:],cx_var[t_start_idx:],label = 'CX Angular Momentum Loss')
axs[1].legend(loc='center left')

power_var = nc_root.variables["Current"][:] * nc_root.variables["Voltage"][:]
ax4 = axs[1].twinx()
ax4.set_ylabel("kW")
ax4.tick_params(axis='y', labelcolor='tab:red')
ax4.plot(t_var[t_start_idx:], power_var[t_start_idx:]/1000,label = 'Power Draw',color = 'tab:red')
ax4.legend(loc='upper right')

# plasmaVolume = np.pi * (0.2 + 2 * 0.15) * 0.2 * 0.6
# multiplier = 2e6
#
# viscous_var = nc_root.variables["ViscousTorque"]
# axs[1].plot(t_var[t_start_idx:],viscous_var[t_start_idx:] * plasmaVolume * multiplier,label = 'Viscous Torque')
# parallel_var = nc_root.variables["ParAngMomLoss"]
# axs[1].plot(t_var[t_start_idx:],parallel_var[t_start_idx:] * plasmaVolume * multiplier,label = 'Parallel Angular Momentum Loss')
# cx_var = nc_root.variables["CXAngMomLoss"]
# axs[1].plot(t_var[t_start_idx:],cx_var[t_start_idx:] * multiplier,label = 'CX Angular Momentum Loss')
# power_var = nc_root.variables["Current"][:] * nc_root.variables["Voltage"][:]
# axs[1].plot(t_var[t_start_idx:],power_var[t_start_idx:],label = 'Power Applied',color = 'tab:red')
# axs[1].legend(loc='upper right')

viscous_var = nc_root.variables["ViscousHeating"]
axs[2].plot(t_var[t_start_idx:],viscous_var[t_start_idx:] / 1000,label = 'Viscous Heating')
parallel_var = nc_root.variables["ParIonHeatLoss"]
axs[2].plot(t_var[t_start_idx:],parallel_var[t_start_idx:] / 1000,label = 'Parallel Ion Heat Loss')
parallel_e_var = nc_root.variables["ParElecHeatLoss"]
axs[2].plot(t_var[t_start_idx:],parallel_e_var[t_start_idx:] / 1000,label = 'Parallel Electron Heat Loss')
cx_var = nc_root.variables["PerpHeatLoss"]
axs[2].plot(t_var[t_start_idx:],cx_var[t_start_idx:] / 1000,label = 'Classical Perpendicular Heat Loss')
axs[2].legend(loc='upper right')
axs[2].set_xlabel("t (ms)")
axs[2].set_ylabel("kW" )
plt.savefig(f'tools/{savefile}.png')


plt.show()
