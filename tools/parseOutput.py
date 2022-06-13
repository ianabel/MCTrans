import pandas as pd
import matplotlib.pyplot as plt
import os

folder = '//wsl$/Ubuntu-20.04/home/nickschw2/programs/MCTrans/parameter_sweep/DT_test_cases'
folder = '/Users/Nick/Downloads/DTTests3'
files = [file for file in os.listdir(folder) if file.endswith('.out')]

columns = ['voltage', 'centralField', 'throatField', 'fieldRatio', 'rotationPower', 'MachAlfven', 'Q_scientific', 'ionTemperature', 'electronTemperature', 'nuStar']
results = pd.DataFrame(columns=columns)

def getResult(filepath):
    result = pd.DataFrame(columns=columns)

    with open(filepath, encoding = 'utf8') as f:
        lines = f.readlines()
        for line in lines:
            line.strip()
            if 'Total potential drop is' in line:
                voltage = float(line.split(' ')[4]) # MV
                units = line.split(' ')[5]
                if 'kV' in units:
                    voltage /= 1000

            elif 'Magnetic Field in the Central Cell' in line:
                centralField = float(line.split(' ')[7]) # T

            elif 'Magnetic Field at the Mirror Throat' in line:
                throatField = float(line.split(' ')[7]) # T

            elif 'Power Required (at the plasma) to support rotation' in line:
                rotationPowerText = line.split(' ')[8].strip() # W
                rotationPower = float(rotationPowerText[:-1])
                units = rotationPowerText[-1]
                if units == 'W':
                    rotationPower /= 1e6
                else:
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

            elif 'ν*' in line:
                nuStar = float(line.split(' ')[3])

        fieldRatio = throatField / centralField
        result = pd.DataFrame([[voltage, centralField, throatField, fieldRatio, rotationPower, MachAlfven, Q_scientific, ionTemperature, electronTemperature, nuStar]], columns=columns)
        return result

for file in files:
    filepath = f'{folder}/{file}'
    result = getResult(filepath)
    results = pd.concat([results, result])

fig, axes = plt.subplots(nrows=6, sharex=True, figsize=(10, 15))

for fieldRatio, B_grp in results.groupby('fieldRatio'):
    B_grp = B_grp.sort_values('voltage')
    B_grp.plot(ax=axes[0], x='voltage', y='rotationPower', logy=True, label=f'R = {fieldRatio:.2f}', legend=False)
    B_grp.plot(ax=axes[1], x='voltage', y='MachAlfven', logy=True, label=f'R = {fieldRatio:.2f}', legend=False)
    B_grp.plot(ax=axes[2], x='voltage', y='Q_scientific', logy=True, label=f'R = {fieldRatio:.2f}', legend=False)
    B_grp.plot(ax=axes[3], x='voltage', y='ionTemperature', logy=True, label=f'R = {fieldRatio:.2f}', legend=False)
    B_grp.plot(ax=axes[4], x='voltage', y='electronTemperature', logy=True, label=f'R = {fieldRatio:.2f}', legend=False)
    B_grp.plot(ax=axes[5], x='voltage', y='nuStar', logy=True, label=f'R = {fieldRatio:.2f}', legend=False)

# Add horizontal line at M_A = 1
axes[1].axhline(y=1, color='r', linestyle='--')

# Add  axis labels
axes[0].set_ylabel('$P_{rot}$ (MW)')
axes[1].set_ylabel('$M_A$')
axes[2].set_ylabel('$Q_{sci}$')
axes[3].set_ylabel('$T_{i}$ (keV)')
axes[4].set_ylabel('$T_{e}$ (keV)')
axes[5].set_ylabel('$\\nu^{*}$ (ions)')

# Add gridlines
axes[0].grid('on')
axes[1].grid('on')
axes[2].grid('on')
axes[3].grid('on')
axes[4].grid('on')
axes[5].grid('on')

# Cut off low Q values
axes[2].set_ylim(bottom=0.1)

plt.xlabel('Voltage $\phi$ (MV)')

# One legend for all plots
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='center right', bbox_to_anchor=(0.92, 0.5))

plt.savefig('results.png', dpi=200)
