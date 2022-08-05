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

t_start_idx = 1

print(t_var[t_start_idx:].size)

fig, (ax1,ax3) = plt.subplots(2,1,sharex=True,figsize=(10,6),dpi=210)

fig.suptitle("Time Traces" )
ax1.set_xlabel("t / seconds")
ax1.set_ylabel("Volts" )
ax1.tick_params(axis='y', labelcolor='tab:blue')

ax1.plot(t_var[t_start_idx:],voltages_var[t_start_idx:],label = 'Applied Voltage')

ax1.legend(loc='upper left')
ax2 = ax1.twinx()

ti_var = nc_root.variables["IonTemperature"];
te_var = nc_root.variables["ElectronTemperature"];


ax2.set_ylabel("keV")
ax2.tick_params(axis='y', labelcolor='tab:blue')
ax2.plot(t_var[t_start_idx:],ti_var[t_start_idx:],label = 'Ion Temperature',color = 'tab:blue')
ax2.plot(t_var[t_start_idx:],te_var[t_start_idx:],label = 'Electron Temperature',color = 'tab:green')
ax2.legend(loc='lower right')

fig.tight_layout()

ax3.set_xlabel("t / seconds")
ax3.set_ylabel("W/m^3" )
ax3.tick_params(axis='y', labelcolor='tab:blue')

viscous_var = nc_root.variables["ViscousTorque"]
#ax3.plot(t_var[t_start_idx:],viscous_var[t_start_idx:],label = 'Viscous Torque')
#parallel_var = nc_root.variables["ParAngMomLoss"]
#ax3.plot(t_var[t_start_idx:125],parallel_var[t_start_idx:125],label = 'Parallel Angular Momentum Loss')
ax3.plot(t_var[t_start_idx:],nc_root.variables["MachNumber"][t_start_idx:],label = 'Mach')
ax3.legend(loc='upper left')

j_var = nc_root.variables["Current"]
ax4 = ax3.twinx()
ax4.set_ylabel("A")
ax4.tick_params(axis='y', labelcolor='tab:blue')
#ax4.plot(t_var[t_start_idx:],j_var[t_start_idx:],label = 'Current',color = 'tab:blue')
ax4.plot(t_var[t_start_idx:],nc_root.variables["AmbipolarPhi"][t_start_idx:],label = 'Phi',color='green')
ax4.legend(loc='lower right')




plt.show()


fig, ax = plt.subplots(1,1,figsize=(10,6),dpi=210)
viscous_var = nc_root.variables["ViscousHeating"]
ax.plot(t_var[t_start_idx:],viscous_var[t_start_idx:],label = 'Viscous Heating')
parallel_var = nc_root.variables["ParIonHeatLoss"]
ax.plot(t_var[t_start_idx:],parallel_var[t_start_idx:],label = 'Parallel Ion Heat Loss')
parallel_e_var = nc_root.variables["ParElecHeatLoss"]
ax.plot(t_var[t_start_idx:],parallel_e_var[t_start_idx:],label = 'Parallel Electron Heat Loss')
cx_var = nc_root.variables["PerpHeatLoss"]
ax.plot(t_var[t_start_idx:],cx_var[t_start_idx:],label = 'Classical Perpendicular Heat Loss')
ax.legend(loc='upper left')


plt.show()

