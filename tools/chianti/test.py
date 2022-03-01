import ChiantiPy.core as ch
import numpy as np
import matplotlib.pyplot as plt

eV2K = 11606
ergsPerSecond2Watts = 1e-7

temperature = np.logspace(1, 4.5, 100) # eV
ne = 1e14 # cm^-3

# Radiative Losses
ionList = ['c_1', 'c_2', 'c_3', 'c_4', 'c_5', 'c_6'] # C VI
rl = ch.radLoss(temperature * eV2K, ne, ionList=ionList)
rate = rl.RadLoss['rate'] * ergsPerSecond2Watts # Watts cm^3

plt.loglog(temperature, rate, label='Carbon ions')
plt.xlabel('Temperature (eV)')
plt.ylabel('Radiative Loss (W cm$^3)$')
plt.legend()
plt.savefig('radLoss.png', dpi=150)

# Ionization Equilibria
ZCarbon = 6
carbonIons = ch.ioneq(ZCarbon)
carbonIons.load()
# Change x-axis to eV
setattr(carbonIons, 'Temperature', carbonIons.Temperature / eV2K)
carbonIons.plot()
plt.xlabel('Temperature (eV)')
plt.savefig('carbonIons.png', dpi=150)

# Hydrogen level populations
hydrogenLevels = ch.ion('h_1', temperature=temperature * eV2K, eDensity=ne)
# Change x-axis to eV
setattr(hydrogenLevels, 'Temperature', hydrogenLevels.Temperature / eV2K)
hydrogenLevels.popPlot()
plt.xlabel('Temperature (eV)')
plt.savefig('hydrogenLevels.png', dpi=150)
