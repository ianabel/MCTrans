#!/usr/bin/python3

from netCDF4 import Dataset
import numpy as np
import sys

filename = ""

if len(sys.argv) == 2:
    filename = sys.argv[1]
else:
    print("Usage: ./GenerateVoltageTrace.py <netcdf filename>")
    sys.exit()

nc_root = Dataset(filename, "w", format="NETCDF4")

TimeDim = nc_root.createDimension("Time",None);
TimeVar = nc_root.createVariable("Time","f8",("Time",));

VoltageVar = nc_root.createVariable("Voltage","f8",("Time",));

# Generate a time trace that goes from t=0 to t=t_1 at V_1, then linearly 
# up to V_2 at t_2, then constant to t_end

t1 = 0.010
t2 = 0.015
t_end = 0.100

V1 = 10000
V2 =100000

N_points = 5001

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
