import pandas as pd
import matplotlib.pyplot as plt
import os

name = 'reactor'
variable = 'electronDensity'
folder = f'/Users/Nick/programs/MCTrans/misc_runs/batch_runs/{name}/vary_{variable}'
# folder = f'/Users/Nick/programs/MCTrans/misc_runs/batch_runs/large_sweep/refined_sweep'
files = [file for file in os.listdir(folder) if file.endswith('.out')]

columns = ['voltage', 'electronDensity', 'centralField', 'throatField', 'fieldRatio', 'rotationPower', 'viscousTorquePower', 'parallelLossPower', 'chargeExchangePower', 'MachAlfven', 'Q_scientific', 'ionTemperature', 'electronTemperature', 'temperatureRatio', 'rho_i', 'a', 'L', 'rmin', 'tau_E', 'Mach', 'rho_star', 'nu_star', 'tripleProduct', 'thermalPower']
results = pd.DataFrame(columns=columns)

# cutoffs
M_A_max = 1.25
rho_star_max = 0.3
nu_star_max = 0.3
thermalPower_min = 10

def getResult(filepath):
    result = pd.DataFrame(columns=columns)

    with open(filepath, encoding = 'utf8') as f:
        voltage, electronDensity, centralField, throatField, fieldRatio, rotationPower, viscousTorquePower, parallelLossPower, chargeExchangePower, MachAlfven, Q_scientific, ionTemperature, electronTemperature, temperatureRatio, rho_i, a, L, rmin, tau_E, Mach, rho_star, nu_star, tripleProduct, thermalPower = None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None
        lines = f.readlines()
        for line in lines:
            line.strip()

            if 'Total potential drop is' in line:
                voltage = float(line.split(' ')[4]) # MV
                units = line.split(' ')[5]
                if 'MV' in units:
                    voltage *= 1000

            elif 'Electron Density is' in line:
                electronDensity = float(line.split(' ')[3]) # m^-3

            elif 'Magnetic Field in the Central Cell' in line:
                centralField = float(line.split(' ')[7]) # T

            elif 'Magnetic Field at the Mirror Throat' in line:
                throatField = float(line.split(' ')[7]) # T

            elif 'Power Required (at the plasma) to support rotation' in line:
                rotationPower = float(line.split(' ')[8])
                units = line.split(' ')[9].strip()
                if units == 'W':
                    rotationPower /= 1e6
                elif units == 'kW':
                    rotationPower /= 1e3
                elif units != 'MW':
                    print(units)

            elif 'Power Loss from viscous torque' in line:
                viscousTorquePower = float(line.split(' ')[7])
                units = line.split(' ')[8].strip()
                if units == 'W':
                    viscousTorquePower /= 1e6
                elif units == 'kW':
                    viscousTorquePower /= 1e3
                elif units != 'MW':
                    print(units)

            elif 'Power Loss from parallel loss' in line:
                parallelLossPower = float(line.split(' ')[8])
                units = line.split(' ')[9].strip()
                if units == 'W':
                    parallelLossPower /= 1e6
                elif units == 'kW':
                    parallelLossPower /= 1e3
                elif units != 'MW':
                    print(units)

            elif 'Power Loss from charge exchange' in line:
                if 'was not included' in line:
                    chargeExchangePower = None
                else:
                    chargeExchangePower = float(line.split(' ')[6])
                    units = line.split(' ')[7].strip()
                    if units == 'W':
                        chargeExchangePower /= 1e6
                    elif units == 'kW':
                        chargeExchangePower /= 1e3
                    elif units != 'MW':
                        print(units)

            elif 'Alfven Mach number is' in line:
                MachAlfven = float(line.split(' ')[4])

            elif 'Q_scientific' in line:
                Q_scientific = float(line.split(' ')[10])

            elif 'Ion Temperature is' in line:
                ionTemperature = float(line.split(' ')[3]) # keV
                units = line.split(' ')[4].strip()
                if units == 'eV':
                    ionTemperature /= 1000
                elif units == 'MeV':
                    ionTemperature *= 1000

            elif 'Electron Temperature is' in line:
                electronTemperature = float(line.split(' ')[3]) # keV
                units = line.split(' ')[4].strip()
                if units == 'eV':
                    electronTemperature /= 1000
                elif units == 'MeV':
                    electronTemperature *= 1000

            elif 'Ion Larmor Radius' in line:
                rho_i = float(line.split(' ')[8])

            elif 'Typical plasma scale lengths' in line:
                a = float(line.split(' ')[6])

            elif 'Energy Confinement Time' in line:
                tau_E = float(line.split(' ')[4])
                units = line.split(' ')[5].strip()
                if units == 'ms':
                    tau_E /= 1000

            elif 'PlasmaLength' in line:
                L = float(line.split(' ')[1])

            elif 'Operating Mach number is' in line:
                Mach = float(line.split(' ')[4])

            elif 'ρ*' in line:
                rho_star = float(line.split(' ')[-1])

            elif 'ν*' in line:
                nu_star = float(line.split(' ')[-3])

            elif 'n T τ' in line:
                tripleProduct = float(line.split(' ')[-4])

            elif 'Total Thermal Power Output is' in line:
                thermalPower = float(line.split(' ')[-2])

            elif 'Inner radius of the plasma' in line:
                rmin = float(line.split(' ')[-2])

        # if rho_i/a < 0.25 and MachAlfven < 1.25 and centralField == 3 and Q_scientific > 3:
        #     print(filepath)

        fieldRatio = throatField / centralField
        temperatureRatio = ionTemperature / electronTemperature
        result = pd.DataFrame([[voltage, electronDensity, centralField, throatField, fieldRatio, rotationPower, viscousTorquePower, parallelLossPower, chargeExchangePower, MachAlfven, Q_scientific, ionTemperature, electronTemperature, temperatureRatio, rho_i, a, L, rmin, tau_E, Mach, rho_star, nu_star, tripleProduct, thermalPower]], columns=columns)
        return result

N = len(files)
for i, file in enumerate(files):
    filepath = f'{folder}/{file}'
    result = getResult(filepath)
    results = pd.concat([results, result])
    print(f'{i+1}/{N}')

# results['rho_star'] = results['rho_i'] / results['a']

# valid_results = results[(results['rho_star'] < rho_star_max) & (results['nu_star'] < nu_star_max) & (results['MachAlfven'] < M_A_max) & (results['thermalPower'] >= thermalPower_min)]
# valid_results.to_csv(f'{folder}/results.csv', columns=columns, index=False)



# Name of parameter and its y-label
if name == 'CMFX':
    plottingParameters = {'rotationPower': {'label': '$P_{rot}$ (MW)',
        'logy': False}, 'MachAlfven': {'label': '$M_A$', 'logy': False},
        'ionTemperature': {'label': '$T_{i}$ (keV)', 'logy': False},
        'temperatureRatio': {'label': '$T_i / T_e$', 'logy': False}}
elif name == 'reactor':
    plottingParameters = {'Q_scientific': {'label': '$Q_{sci}$', 'logy': False},
        'ionTemperature': {'label': '$T_{i}$ (keV)', 'logy': False},
        'temperatureRatio': {'label': '$T_i / T_e$', 'logy': False},
        'thermalPower': {'label': '$P_{thermal}$ (MW)', 'logy': True}}
elif name == 'physics':
    plottingParameters = {'tau_E': {'label': r'$\tau_E$', 'logy': False},
        'rho_star': {'label': r'$\rho^*$', 'logy': False}}

fig, axes = plt.subplots(nrows=len(plottingParameters), sharex=True, figsize=(10, 15))
# Change the font size for all plots
plt.rcParams.update({'font.size': 18})

for value, grp in results.groupby(variable):
    grp = grp.sort_values('voltage')

    if name == 'reactor':
        grp['voltage'] = grp['voltage'] / 1000

    # Filter results for plotting
    grp = grp[(grp['rho_star'] < rho_star_max) & (grp['MachAlfven'] < M_A_max)]
    # grp = grp[(grp['rho_star'] < rho_star_max * 1.2)]

    for i, parameter in enumerate(plottingParameters):
        if variable == 'centralField':
            label = f'$B_{{mid}}$ = {value:.2f} T'
        elif variable == 'electronDensity':
            label = f'$n_{{e}}$ = {value:.0e} m$^{{-3}}$'

        grp.plot(ax=axes[i], x='voltage', y=parameter, logy=plottingParameters[parameter]['logy'], label=label, legend=False, fontsize=12)
        ylabel = plottingParameters[parameter]['label']
        axes[i].set_ylabel(ylabel, fontsize=16)
        axes[i].grid('on')
        axes[i].tick_params(axis='both', labelsize=16)

        if parameter == 'MachAlfven':
            axes[i].axhline(y=1, color='r', linestyle='--')
        elif parameter == 'rho_star':
            axes[i].axhline(y=rho_star_max, color='r', linestyle='--')
        elif parameter == 'thermalPower':
            # Cut off low thermal power values
            axes[i].set_ylim(bottom=0.1, auto=True)

if name == 'reactor':
    plt.xlabel('Voltage $\phi$ (MV)', fontsize=16)
else:
    plt.xlabel('Voltage $\phi$ (kV)', fontsize=16)

# One legend for all plots
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='center left', bbox_to_anchor=(0.13, 0.81), prop={'size': 12})
plt.savefig(f'results_{name}_{variable}.png', dpi=400)
plt.show()

# Plot different sources of power draw for one curve
grp = results[results[variable] == min(results[variable])]
grp = grp.sort_values('voltage')

if name == 'reactor':
    grp['voltage'] = grp['voltage'] / 1000

parameters = {'rotationPower': 'Total', 'viscousTorquePower': 'Viscous Torque', 'parallelLossPower': '|| Loss', 'chargeExchangePower': 'Charge Exchange'}

fig, axes = plt.subplots(figsize=(8, 6))
axes.tick_params(axis='both', labelsize=16)
for key in parameters:
    label=parameters[key]
    if name == 'CMFX':
        grp[key] = grp[key] * 1000
    grp.plot(ax=axes, x='voltage', y=key, label=label)

if name == 'reactor':
    plt.xlabel('Voltage $\phi$ (MV)', fontsize=16)
    plt.ylabel('$P$ (MW)')
else:
    plt.xlabel('Voltage $\phi$ (kV)', fontsize=16)
    plt.ylabel('$P$ (kW)')
plt.legend()
plt.savefig(f'results_power_{name}_{variable}.png', dpi=400)
plt.show()

# Graph for physics understanding
physicsFolder = '/Users/Nick/programs/MCTrans/misc_runs/batch_runs/physics'
physicsFiles = [file for file in os.listdir(physicsFolder) if file.endswith('.out')]
physicsResults = pd.DataFrame(columns=columns)

CXBool = [True] * len(physicsFiles)
PhiBool = [True] * len(physicsFiles)
for i, file in enumerate(physicsFiles):
    if 'noCX' in file:
        CXBool[i] = False
    if 'noPhi' in file:
        PhiBool[i] = False

    filepath = f'{physicsFolder}/{file}'
    result = getResult(filepath)
    physicsResults = pd.concat([physicsResults, result])

# Add columns to dataframe
physicsResults['CX'] = CXBool
physicsResults['Phi'] = PhiBool

physicsResults = physicsResults.sort_values('voltage')
resultsNormal = physicsResults[(physicsResults['Phi'] == True) & (physicsResults['CX'] == True)]
resultsNoCX = physicsResults[(physicsResults['CX'] == False)]
resultsNoPhi = physicsResults[(physicsResults['Phi'] == False)]

plottingParameters = {'rotationPower': '$P_{rot}$ (MW)', 'ionTemperature': '$T_{i}$ (keV)', 'temperatureRatio': '$T_i / T_e$'}
fig, axes = plt.subplots(nrows=len(plottingParameters), sharex=True, figsize=(10, 15))

for i, parameter in enumerate(plottingParameters):
    resultsNormal.plot(ax=axes[i], x='voltage', y=parameter, logy=False, label='Normal', legend=False, fontsize=12)
    resultsNoCX.plot(ax=axes[i], x='voltage', y=parameter, logy=False, label='No CX', legend=False, fontsize=12)
    resultsNoPhi.plot(ax=axes[i], x='voltage', y=parameter, logy=False, label=r'No $\varphi$', legend=False, fontsize=12)
    ylabel = plottingParameters[parameter]
    axes[i].set_ylabel(ylabel, fontsize=16)
    axes[i].grid('on')
    # axes[i].tick_params(axis='both', labelsize=16)

    if parameter == 'MachAlfven':
        axes[i].axhline(y=1, color='r', linestyle='--')
    elif parameter == 'Q_scientific':
        # Cut off low Q values
        axes[i].set_ylim(bottom=0.1)

plt.xlabel('Voltage $\phi$ (kV)', fontsize=16)

# One legend for all plots
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='center left', bbox_to_anchor=(0.13, 0.8))
plt.savefig(f'physics_results.png', dpi=400)
plt.show()
