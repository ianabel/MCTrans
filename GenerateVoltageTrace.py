#!/usr/bin/env python3

from netCDF4 import Dataset
import numpy as np
import sys

filename = ""

# Generate a time trace that goes from t=0 to t=t_1 at V_1, then linearly 
# up to V_2 at t_2, then constant to t_end

if len(sys.argv) == 7:
    filename = sys.argv[1]
    V1 = float(sys.argv[2])
    V2 = float(sys.argv[3])
    t1 = float(sys.argv[4])
    t2 = float(sys.argv[5])
    t_end = float(sys.argv[6])
else:
    print("Usage: ./GenerateVoltageTrace.py <netcdf filename> <start voltage> <end voltage> <ramp start> <ramp end> <end time>")
    sys.exit()

nc_root = Dataset(filename, "w", format="NETCDF4")

TimeDim = nc_root.createDimension("Time",None);
TimeVar = nc_root.createVariable("Time","f8",("Time",));

VoltageVar = nc_root.createVariable("Voltage","f8",("Time",));


N_points = 15001

delta_t = t_end / ( N_points - 1 )

for i in range(0,N_points):
    time = delta_t * i;
    TimeVar[i] = time;
    if time <= t1:
        VoltageVar[i] = V1;
    if time > t1 and time <= t2:
        VoltageVar[i] = V1 + ( V2 - V1 ) * ( time - t1 ) / ( t2 - t1 )
    if time > t2:
        VoltageVar[i] = V2

nc_root.close()
