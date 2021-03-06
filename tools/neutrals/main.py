import matplotlib.pyplot as plt
from scipy.integrate import cumtrapz
from neutrals import *

temperature = np.logspace(1, 5, 100) # eV
energy = np.logspace(1, 6, 200) # eV
density = 1e20 # m^-3
plasmaVolume = 10 # m^3
MachNumber = 4

neutrals = Neutrals(MachNumber, temperature, energy)

# Plot cross sections
plt.figure(1)
for crossSection in neutrals.crossSections :
    plt.loglog(energy, crossSection.sigma, label=crossSection.reaction)

plt.xlabel('Energy (eV)')
plt.ylabel('$\sigma_{cx}$ (cm$^2$)')
plt.legend()
plt.grid()
plt.savefig('CrossSections.png', dpi=150)

# Plot rate coefficients
prop_cycle = plt.rcParams['axes.prop_cycle']
colors = prop_cycle.by_key()['color']
for color, rateCoefficientHot, rateCoefficientCold in zip(colors, neutrals.rateCoefficientsHot, neutrals.rateCoefficientsCold):
    plt.figure(2)
    plt.loglog(temperature, rateCoefficientHot.rateCoefficient, label=rateCoefficientHot.reaction)

    plt.figure(3)
    plt.loglog(temperature, rateCoefficientHot.rateCoefficient.T, color=color, label=rateCoefficientCold.reaction)
    plt.loglog(temperature, rateCoefficientCold.rateCoefficient.T, '--', color=color)

plt.figure(2)
plt.title('Hot Rate Coefficients')
plt.xlabel('T (eV)')
plt.ylabel('$<\sigma v>$ (cm$^3$/s)')
plt.legend()
plt.grid()
plt.ylim([1e-14, 1e-6])
plt.xlim([1e1, 1e5])
plt.savefig('HotRateCoefficients.png', dpi=150)

plt.figure(3)
plt.title(f'Cold Rate Coefficients (M={MachNumber})')
plt.xlabel('T (eV)')
plt.ylabel('$<\sigma v>$ (cm$^3$/s)')
plt.legend()
plt.grid()
plt.ylim([1e-14, 1e-6])
plt.xlim([1e1, 1e5])
plt.savefig('ColdRateCoefficients.png', dpi=150)
