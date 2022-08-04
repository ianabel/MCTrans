#!/usr/bin/env python3

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

t1 = 0.0050
t2 = 0.0550
t_end = 0.0750

V1 = 10000
V2 =100000

N_points = 15001

delta_t = t_end / ( N_points - 1 )

def sigmoid(x):
    sig = np.where(x < 0, np.exp(x)/(1 + np.exp(x)), 1/(1 + np.exp(-x)))
    return sig


for i in range(0,N_points):
    time = delta_t * i;
    TimeVar[i] = time;
    VoltageVar[i] = V1 + ( V2 - V1 ) * sigmoid( 10*(time - (t1+t2)/2.0)/(t2-t1) );

nc_root.close()