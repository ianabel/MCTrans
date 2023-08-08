import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import numpy as np
from concurrent.futures import ProcessPoolExecutor

# Change color cycler
# cmap = plt.get_cmap('plasma')
# color_palette = cmap(np.linspace(0, 1, 7))[0:6]
# color_palette = ['#15005c', '#6e0061', '#ae0059', '#dd2646', '#f8692c', '#ffa600']
# mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=color_palette)
# Change the font size for all plots
mpl.rcParams.update({'lines.linewidth': 2})
labelsize = 15
reset = False

name = 'reactor'
variable = 'electronDensity'
folder = f'/Users/Nick/programs/MCTrans/misc_runs/batch_runs/{name}/vary_{variable}'
folder = f'/Users/Nick/programs/MCTrans/misc_runs/batch_runs/large_sweep'
files = [f'{folder}/{file}' for file in os.listdir(folder) if file.endswith('.out')]

# for filename in files:
#     nc_root = Dataset(f'{folder}/{filename}', "r", format="NETCDF4")
#     breakpoint()

columns = ['voltage',
           'electronDensity',
           'centralField',
           'throatField',
           'fieldRatio',
           'rotationPower',
           'viscousTorquePower',
           'parallelLossPower',
           'chargeExchangePower',
           'totalHeatLoss',
           'classicalIonHeatLoss',
           'parallelIonHeatLoss',
           'chargeExchangeHeatLoss',
           'parallelElectronHeatLoss',
           'radiationLoss',
           'MachAlfven',
           'Q_scientific',
           'Q_scientific_equiv',
           'ionTemperature',
           'electronTemperature',
           'temperatureRatio',
           'rho_i',
           'a',
           'L',
           'rmin',
           'tau_E',
           'Mach',
           'rho_star',
           'nu_star',
           'tripleProduct',
           'thermalPower']

# cutoffs
M_A_max = 1.25
rho_star_max = 0.3
nu_star_max = 0.3
thermalPower_min = 50 # MW
Qsci_min = 1
voltage_max = 1000

def getResult(filepath):
    result = pd.DataFrame(columns=columns)

    with open(filepath, encoding = 'utf8') as f:
        voltage, electronDensity, centralField, throatField, fieldRatio, rotationPower, viscousTorquePower, parallelLossPower, chargeExchangePower, totalHeatLoss, classicalIonHeatLoss, parallelIonHeatLoss, chargeExchangeHeatLoss, parallelElectronHeatLoss, radiationLoss, MachAlfven, Q_scientific, Q_scientific_equiv, ionTemperature, electronTemperature, temperatureRatio, rho_i, a, L, rmin, tau_E, Mach, rho_star, nu_star, tripleProduct, thermalPower = None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None
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
                    rotationPower /= 1e3
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

            elif 'Total heat losses is' in line:
                totalHeatLoss = float(line.split(' ')[4])
                units = line.split(' ')[5].strip()
                if units == 'W':
                    totalHeatLoss /= 1e6
                elif units == 'kW':
                    totalHeatLoss /= 1e3
                elif units != 'MW':
                    print(units)

            elif 'Heat Loss from classical ion loss' in line:
                classicalIonHeatLoss = float(line.split(' ')[8])
                units = line.split(' ')[9].strip()
                if units == 'W':
                    classicalIonHeatLoss /= 1e6
                elif units == 'kW':
                    classicalIonHeatLoss /= 1e3
                elif units != 'MW':
                    print(units)

            elif 'Heat Loss from parallel ion loss' in line:
                parallelIonHeatLoss = float(line.split(' ')[8])
                units = line.split(' ')[9].strip()
                if units == 'W':
                    parallelIonHeatLoss /= 1e6
                elif units == 'kW':
                    parallelIonHeatLoss /= 1e3
                elif units != 'MW':
                    print(units)

            elif 'Heat loss due to charge exchange' in line:
                if 'was not included' in line:
                    chargeExchangeHeatLoss = None
                else:
                    chargeExchangeHeatLoss = float(line.split(' ')[10])
                    units = line.split(' ')[11].strip()
                    if units == 'W':
                        chargeExchangeHeatLoss /= 1e6
                    elif units == 'kW':
                        chargeExchangeHeatLoss /= 1e3
                    elif units != 'MW':
                        print(units)

            elif 'Heat Loss from parallel electron loss' in line:
                parallelElectronHeatLoss = float(line.split(' ')[8])
                units = line.split(' ')[9].strip()
                if units == 'W':
                    parallelElectronHeatLoss /= 1e6
                elif units == 'kW':
                    parallelElectronHeatLoss /= 1e3
                elif units != 'MW':
                    print(units)

            elif 'Heat Loss from radiation' in line:
                radiationLoss = float(line.split(' ')[6])
                units = line.split(' ')[7].strip()
                if units == 'W':
                    radiationLoss /= 1e6
                elif units == 'kW':
                    radiationLoss /= 1e3
                elif units != 'MW':
                    print(units)

            elif 'Alfven Mach number is' in line:
                MachAlfven = float(line.split(' ')[4])

            elif 'Q_scientific' in line:
                Q_scientific = float(line.split(' ')[10])

            elif ' Q_(Scientific DT equivalent)' in line:
                Q_scientific_equiv = float(line.split(' ')[5])

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
        result = pd.DataFrame([[voltage, electronDensity, centralField, throatField, fieldRatio, rotationPower, viscousTorquePower, parallelLossPower, chargeExchangePower, totalHeatLoss, classicalIonHeatLoss, parallelIonHeatLoss, chargeExchangeHeatLoss, parallelElectronHeatLoss, radiationLoss, MachAlfven, Q_scientific, Q_scientific_equiv, ionTemperature, electronTemperature, temperatureRatio, rho_i, a, L, rmin, tau_E, Mach, rho_star, nu_star, tripleProduct, thermalPower]], columns=columns)
        return result

if __name__ == '__main__':
    results = pd.DataFrame(columns=columns)
    results_list = []
    if reset:
        N = len(files)
        with ProcessPoolExecutor(max_workers=4) as executor:
            for i, result in enumerate(executor.map(getResult, files)):
                results_list.append(result)
                print(f'{i+1}/{N}')
        results = pd.concat(results_list)
        results['rho_star'] = results['rho_i'] / results['a']

        results.to_csv(f'{folder}/results.csv', columns=columns, index=False)

    else:
        results = pd.read_csv(f'{folder}/results.csv')

    valid_results = results[(results['rho_star'] < rho_star_max) & (results['MachAlfven'] < M_A_max) & (results['voltage'] <= voltage_max)]
    valid_results.to_csv(f'{folder}/results_filtered.csv', columns=columns, index=False)

    # # Name of parameter and its y-label
    # if name == 'CMFX':
    #     plottingParameters = {'rotationPower': {'label': '$P_{in}$ (kW)','logy': False},
    #         'MachAlfven': {'label': '$M_A$', 'logy': False},
    #         'ionTemperature': {'label': '$T_i$ (keV)', 'logy': False},
    #         'temperatureRatio': {'label': '$T_i / T_e$', 'logy': False}}
    # elif name == 'reactor':
    #     plottingParameters = {'Q_scientific': {'label': '$Q_{sci}$', 'logy': False},
    #         'thermalPower': {'label': '$P_{thermal}$ (MW)', 'logy': True},
    #         'ionTemperature': {'label': '$T_{i}$ (keV)', 'logy': False},
    #         'temperatureRatio': {'label': '$T_i / T_e$', 'logy': False}}
    # elif name == 'physics':
    #     plottingParameters = {'tau_E': {'label': r'$\tau_E$', 'logy': False},
    #         'rho_star': {'label': r'$\rho^*$', 'logy': False}}

    # mpl.rcParams.update({'font.size': 20})
    # fig, axes = plt.subplots(nrows=len(plottingParameters), sharex=True, figsize=(8, 16))

    # # Initialize dataframes to show limiting lines
    # M_A_points = pd.DataFrame(columns=columns)
    # rho_star_points = pd.DataFrame(columns=columns)

    # # We don't want to plot every single line, so only plot up to n_lines
    # n_lines = 6
    # n_variable = len(results.groupby(variable))
    # i_list = np.linspace(0, n_variable - n_variable % n_lines, n_lines).astype('int')
    # if variable == 'electronDensity':
    #     if name == 'CMFX':
    #         i_list = [0, 1, 4, 9, 13, 18]
    #     elif name == 'reactor':
    #         i_list = [0, 4, 9, 13, 18, 28]
    # elif variable == 'centralField':
    #     if name == 'CMFX':
    #         i_list = [0, 5, 10, 15, 20]
    #     elif name == 'reactor':
    #         i_list = [10, 13, 16, 19, 22, 25]

    # for i, (value, grp) in enumerate(results.groupby(variable)):
    #     grp = grp.sort_values('voltage')

    #     if name == 'reactor':
    #         grp['voltage'] = grp['voltage'] / 1000

    #     if name == 'CMFX':
    #         grp['rotationPower'] = grp['rotationPower'] * 1000

    #     # Add lines where limits are crossed by finding first row at which limit is crossed
    #     # MachAlfven first
    #     M_A_filter = grp[grp['MachAlfven'] > M_A_max]
    #     rho_star_filter = grp[grp['rho_star'] > rho_star_max]
    #     if len(M_A_filter) > 0:
    #         M_A_point = M_A_filter.iloc[0]
    #         M_A_points = pd.concat([M_A_points, M_A_point.to_frame().T])
    #     if len(rho_star_filter) > 0:
    #         rho_star_point = rho_star_filter.iloc[0]
    #         rho_star_points = pd.concat([rho_star_points, rho_star_point.to_frame().T])

    #     # Only plot every nth loop...
    #     if i not in i_list:
    #         continue

    #     # Filter results for plotting
    #     grp = grp[(grp['rho_star'] < rho_star_max) & (grp['MachAlfven'] < M_A_max)]

    #     for i, parameter in enumerate(plottingParameters):
    #         if variable == 'centralField':
    #             label = f'$B_{{min}}$ = {value:.2f} T'
    #         elif variable == 'electronDensity':
    #             label = f'$n_{{e}}$ = {value:.0e} m$^{{-3}}$'

    #         grp.plot(ax=axes[i], x='voltage', y=parameter, logy=plottingParameters[parameter]['logy'], label=label, legend=False)

    # # Plot the limiting lines
    # for i, parameter in enumerate(plottingParameters):
    #     if variable == 'electronDensity':
    #         rho_star_points.plot(style='g--', ax=axes[i], x='voltage', y=parameter, label = '_nolegend_', logy=plottingParameters[parameter]['logy'], legend=False)
    #         if parameter == 'MachAlfven':
    #             axes[i].axhline(y=M_A_max, color='m', linestyle='--')
    #         elif parameter != 'ionTemperature' and parameter != 'temperatureRatio' and parameter != 'Q_scientific':
    #             M_A_points.plot(style='m--', ax=axes[i], x='voltage', y=parameter, label = '_nolegend_', logy=plottingParameters[parameter]['logy'], legend=False)
    #     elif variable == 'centralField':
    #         rho_star_points.plot(style='g--', ax=axes[i], x='voltage', y=parameter, label = '_nolegend_', logy=plottingParameters[parameter]['logy'], legend=False)
    #         if parameter == 'MachAlfven':
    #             axes[i].axhline(y=M_A_max, color='m', linestyle='--')
    #         else:
    #             M_A_points.plot(style='m--', ax=axes[i], x='voltage', y=parameter, label = '_nolegend_', logy=plottingParameters[parameter]['logy'], legend=False)
        
    #     ylabel = plottingParameters[parameter]['label']
    #     axes[i].set_ylabel(ylabel)
    #     axes[i].grid('on')
    #     axes[i].tick_params(axis='both', labelsize=labelsize)

    #     # if parameter == 'MachAlfven':
    #     #     axes[i].axhline(y=1, color='r', linestyle='--')
    #     # elif parameter == 'rho_star':
    #     #     axes[i].axhline(y=rho_star_max, color='r', linestyle='--')
    #     if parameter == 'thermalPower':
    #         # Cut off low thermal power values
    #         axes[i].set_ylim(bottom=1.0, auto=True)

    # if name == 'reactor':
    #     plt.xlabel('Voltage $\phi$ (MV)')
    #     # if variable == 'electronDensity':
    #         # plt.xlim([0.8, 10])
    #         # axes[2].set_ylim(top=180)
    #         # axes[3].set_ylim(top=10)
    # elif name == 'CMFX':
    #     plt.xlabel('Voltage $\phi$ (kV)')
    #     if variable == 'electronDensity':
    #         axes[0].set_ylim(top=1950)
    #         axes[1].set_ylim(top=1.3)

    # # One legend for all plots
    # handles, labels = axes[0].get_legend_handles_labels()
    # bbox_to_anchor = (0.13, 0.81)
    # if name == 'CMFX':
    #     if variable == 'electronDensity':
    #         bbox_to_anchor = (0.62, 0.81)
    #     elif variable == 'centralField':
    #         bbox_to_anchor = (0.65, 0.83)
    # fig.legend(handles, labels, loc='center left', bbox_to_anchor=bbox_to_anchor, prop={'size': 12})
    # plt.savefig(f'results_{name}_{variable}.png', dpi=400)
    # plt.show()

    # ### POWER COMPARISON ###
    # folder = f'/Users/Nick/programs/MCTrans/misc_runs/batch_runs/power_comparison/{name}'
    # files = [file for file in os.listdir(folder) if file.endswith('.out')]
    # results = pd.DataFrame(columns=columns)
    # N = len(files)
    # for i, file in enumerate(files):
    #     filepath = f'{folder}/{file}'
    #     result = getResult(filepath)
    #     results = pd.concat([results, result])
    #     print(f'{i+1}/{N}')

    # results['rho_star'] = results['rho_i'] / results['a']

    # results.to_csv(f'{folder}/results.csv', columns=columns, index=False)

    # # Don't plot non-physical results
    # results = results[(results['MachAlfven'] <= M_A_max) & (results['rho_star'] <= rho_star_max)]

    # # Plot different sources of power draw for one curve
    # results = results.sort_values('voltage')

    # if name == 'reactor':
    #     results['voltage'] = results['voltage'] / 1000

    # heat_parameters = {'totalHeatLoss': 'Total', 'classicalIonHeatLoss': '$\perp$ Ion Loss', 'parallelIonHeatLoss': '|| Ion Loss', 'chargeExchangeHeatLoss': 'Charge Exchange', 'parallelElectronHeatLoss': '|| Electron Loss', 'radiationLoss': 'Radiation'}
    # momentum_parameters = {'rotationPower': 'Total', 'viscousTorquePower': 'Viscous Torque', 'parallelLossPower': '|| Ion Loss', 'chargeExchangePower': 'Charge Exchange'}

    # mpl.rcParams.update({'font.size': 15})
    # fig, axes = plt.subplots(nrows=2, sharex=True, figsize=(8, 8))
    # for key, label in momentum_parameters.items():
    #     if name == 'CMFX':
    #         results[key] = results[key] * 1000
    #     results.plot(ax=axes[0], x='voltage', y=key, label=label)

    # for key, label in heat_parameters.items():
    #     if name == 'CMFX':
    #         results[key] = results[key] * 1000
    #     results.plot(ax=axes[1], x='voltage', y=key, label=label)

    # if name == 'reactor':
    #     plt.xlabel('Voltage $\phi$ (MV)')
    #     axes[0].set_ylabel('Momentum Loss (MW)')
    #     axes[1].set_ylabel('Heat Loss (MW)')
    # else:
    #     plt.xlabel('Voltage $\phi$ (kV)')
    #     axes[0].set_ylabel('Momentum Loss (kW)')
    #     axes[1].set_ylabel('Heat Loss (kW)')
    # plt.legend()
    # axes[0].grid('on')
    # axes[1].grid('on')
    # axes[0].tick_params(axis='both', labelsize=labelsize)
    # axes[1].tick_params(axis='both', labelsize=labelsize)
    # plt.savefig(f'results_losses_{name}.png', dpi=400)
    # plt.show()

    # ### PHYSICS UNDERSTANDING ###
    # physicsFolder = '/Users/Nick/programs/MCTrans/misc_runs/batch_runs/physics'
    # physicsFiles = [file for file in os.listdir(f'{physicsFolder}/{name}') if file.endswith('.out')]
    # physicsResults = pd.DataFrame(columns=columns)

    # CXBool = [True] * len(physicsFiles)
    # PhiBool = [True] * len(physicsFiles)
    # for i, file in enumerate(physicsFiles):
    #     if 'noCX' in file:
    #         CXBool[i] = False
    #     if 'noPhi' in file:
    #         PhiBool[i] = False

    #     filepath = f'{physicsFolder}/{name}/{file}'
    #     result = getResult(filepath)
    #     physicsResults = pd.concat([physicsResults, result])

    # # Add columns to dataframe
    # physicsResults['CX'] = CXBool
    # physicsResults['Phi'] = PhiBool

    # if name == 'reactor':
    #     physicsResults['voltage'] = physicsResults['voltage'] / 1000

    # physicsResults = physicsResults.sort_values('voltage')
    # # Cut out non-physical results
    # physicsResults = physicsResults[(physicsResults['MachAlfven'] <= M_A_max) & (physicsResults['rho_star'] <= rho_star_max)]
    # physicsResults['rotationPower'] *= 1000
    # resultsNormal = physicsResults[(physicsResults['Phi'] == True) & (physicsResults['CX'] == True)]
    # resultsNoCX = physicsResults[(physicsResults['CX'] == False)]
    # resultsNoPhi = physicsResults[(physicsResults['Phi'] == False)]

    # mpl.rcParams.update({'font.size': 20})
    # fig, axes = plt.subplots(nrows=len(plottingParameters), sharex=True, figsize=(8, 16))

    # for i, parameter in enumerate(plottingParameters):
    #     resultsNormal.plot(ax=axes[i], x='voltage', y=parameter, logy=False, label='Normal', legend=False)
    #     resultsNoCX.plot(ax=axes[i], x='voltage', y=parameter, logy=False, label='No CX', legend=False)
    #     resultsNoPhi.plot(ax=axes[i], x='voltage', y=parameter, logy=False, label=r'No $\varphi$', legend=False)
    #     ylabel = plottingParameters[parameter]['label']
    #     axes[i].set_ylabel(ylabel)
    #     axes[i].grid('on')
    #     axes[i].tick_params(axis='both', labelsize=labelsize)

    #     # if parameter == 'MachAlfven':
    #     #     axes[i].axhline(y=1, color='r', linestyle='--')
    #     # elif parameter == 'Q_scientific':
    #     #     # Cut off low Q values
    #     #     axes[i].set_ylim(bottom=0.1)

    # if name == 'reactor':
    #     plt.xlabel('Voltage $\phi$ (MV)')
    # else:
    #     plt.xlabel('Voltage $\phi$ (kV)')

    # # One legend for all plots
    # handles, labels = axes[0].get_legend_handles_labels()
    # fig.legend(handles, labels, loc='center left', bbox_to_anchor=(0.13, 0.8), prop={'size': 16})
    # plt.savefig(f'physics_results_{name}.png', dpi=400)
    # plt.show()
