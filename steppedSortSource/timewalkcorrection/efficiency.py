import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Data from the provided table
distance = np.array([10]*9 + [30]*9 + [50]*9 + [100]*9 + [200]*9)
solidangle = np.array([2.5183]*9 + [1.4871]*9 + [0.8111]*9 + [0.3238]*9 + [0.0881]*9)

energy = np.array([
    121.78, 244.69, 344.28, 411.12, 443.98, 778.9, 867.4, 964.13, 1408,
    121.78, 244.69, 344.28, 411.12, 443.98, 778.9, 867.4, 964.13, 1408,
    121.78, 244.69, 344.28, 411.12, 443.98, 778.9, 867.4, 964.13, 1408,
    121.78, 244.69, 344.28, 411.12, 443.98, 778.9, 867.4, 964.13, 1408,
    121.78, 244.69, 344.28, 411.12, 443.98, 778.9, 867.4, 964.13, 1408
])
fepe = np.array([
    0.0630, 0.0457, 0.0399, 0.0332, 0.0308, 0.0170, 0.0117, 0.0090, 0.0062,
    0.0407, 0.0296, 0.0248, 0.0214, 0.0190, 0.0104, 0.0089, 0.0080, 0.0071,
    0.0312, 0.0202, 0.0177, 0.0147, 0.0132, 0.0075, 0.0065, 0.0060, 0.0059,
    0.0102, 0.0086, 0.0070, 0.0058, 0.0055, 0.0038, 0.0031, 0.0026, 0.0020,
    0.0029, 0.0025, 0.0021, 0.0019, 0.0017, 0.0011, 0.0009, 0.0008, 0.0005
])
log_energy = np.log10(energy)
err_fepe = np.array([
    0.0071, 0.0051, 0.0045, 0.0037, 0.0035, 0.0020, 0.0014, 0.0011, 0.0007,
    0.0044, 0.0030, 0.0025, 0.0022, 0.0019, 0.0011, 0.0009, 0.0008, 0.0007,
    0.0042, 0.0020, 0.0018, 0.0015, 0.0013, 0.0008, 0.0007, 0.0006, 0.0006,
    0.0021, 0.0009, 0.0007, 0.0008, 0.0006, 0.0004, 0.0003, 0.0002, 0.0003,
    0.0004, 0.0004, 0.0006, 0.0006, 0.0004, 0.0002, 0.0003, 0.0002, 0.0001
])
fepeint = np.array([
    0.31409493, 0.22776497, 0.19861872, 0.16535384, 0.15323379, 0.08448635,
0.05805218, 0.04483254, 0.03091594, 0.34849103, 0.25373904,
0.21255136, 0.18360389, 0.16256884, 0.08905580, 0.07647527,
0.06889618, 0.06120046, 0.37037484, 0.30137576, 0.26435544, 0.21853394,
0.19645560, 0.11172307, 0.09695785, 0.09011255, 0.08822440, 0.38477491,
0.32302102, 0.26375788, 0.21943300, 0.20779955, 0.14177385, 0.11632442,
0.09821390, 0.07680860, 0.39363128, 0.34498813, 0.29047328, 0.25754718,
0.23781329, 0.14721639, 0.12098225, 0.11193478, 0.07083656
])

# Unique distances
distances = np.unique(distance)
colors = ['b', 'g', 'r', 'c', 'm']

# Make axis labels bold
plt.rcParams["axes.labelweight"] = "bold"
# Make all plots square
#make all fonts 22
plt.rcParams.update({'font.size': 28})
plt.tick_params(axis='both', which='major', labelsize=26)
plt.tick_params(axis='both', which='major', labelsize=26)

plt.rcParams['figure.figsize'] = [10, 10]

def log_fit(x, a, b):
    return a * np.log(x) + b

def linear_fit(x, a, b):
    return a * x + b    

# Define the double exponential decay function
def double_exp_decay(x, a, b, c, d):
    return a * np.exp(-x / b) + c * np.exp(-x / d)

# Plot settings for a publication
plt.rcParams.update({'font.size': 12, 'figure.figsize': (8, 6), 'lines.linewidth': 1.5, 'grid.color': 'grey', 'grid.linestyle': '--', 'grid.linewidth': 0.5})

# Plot individual graphs of Energy vs FEPE and FEPEint with blue circle markers
# Plot individual graphs of Energy vs FEPE for each distance
""" for i, dist in enumerate(distances):
    mask = distance == dist
    x = energy[mask]
    y = fepe[mask]

    plt.figure()
    plt.errorbar(x, y, yerr=err_fepe[mask], fmt='o', color='blue', ecolor='blue', elinewidth=0.75, capsize=3, label=f'{dist} mm')
    popt, _ = curve_fit(log_fit, x, y)
    linspace = np.linspace(min(x), max(x), 100)
    plt.plot(linspace, log_fit(linspace, *popt), color='blue', linestyle='--', linewidth=1)
    plt.xlabel('Energy (keV)')
    plt.ylabel(r'$\mathbf{\epsilon}_{\mathbf{FEP}}$')
    plt.legend(loc='upper right')
    plt.grid(True)
    plt.tight_layout()    
    plt.show()  """

# Plot individual graphs of Energy vs FEPEint for each distance
""" for i, dist in enumerate(distances):
    mask = distance == dist
    x = energy[mask]
    yint = fepeint[mask]

    plt.figure()
    plt.errorbar(x, yint, yerr=err_fepe[mask], fmt='o', color='blue', ecolor='blue', elinewidth=0.75, capsize=3, label=f'{dist} mm')
    popt, _ = curve_fit(log_fit, x, yint)
    linspace = np.linspace(min(x), max(x), 100)
    plt.plot(linspace, log_fit(linspace, *popt), color='blue', linestyle='--', linewidth=1)
    plt.xlabel('Energy (keV)')
    plt.ylabel(r'$\mathbf{\epsilon}_{\mathbf{FEP,int}}$')    
    plt.legend(loc='upper right')
    plt.grid(True)
    plt.tight_layout()    
    plt.show()  """


markers = ['o', 's', 'D', '^', 'v']
colours = ['b', 'g', 'r', 'c', 'm']

""" # Plot of FEPE vs log(Energy) for all distances
plt.figure()
for i, dist in enumerate(distances):
    mask = distance == dist
    x = log_energy[mask]
    y = fepe[mask]

    plt.plot(x, y, markers[i], color=colours[i],label=f'{dist} mm')
    #linear fit
    popt, _ = curve_fit(linear_fit, x, y)
    plt.plot(x, linear_fit(x, *popt), color=colours[i], linestyle='--', linewidth=1)
plt.xlabel('log(Energy) (log(keV))')
plt.ylabel(r'$\mathbf{\epsilon}_{\mathbf{FEP}}$')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()  """ 

""" # Plot of FEPEint vs log(Energy) for all distances
plt.figure()
for i, dist in enumerate(distances):
    mask = distance == dist
    x = log_energy[mask]
    yint = fepeint[mask]

    plt.plot(x, yint,  markers[i], color=colours[i],label=f'{dist} mm')
    #linear fit
    popt, _ = curve_fit(linear_fit, x, yint)
    plt.plot(x, linear_fit(x, *popt), color=colours[i], linestyle='--', linewidth=1)

plt.xlabel('log(Energy) (log(keV))')
plt.ylabel(r'$\mathbf{\epsilon}_{\mathbf{FEP,int}}$')
plt.grid(True)
plt.ylim([0, 0.41])
plt.legend()
plt.tight_layout()
plt.show() """

#make all fonts 28
plt.rcParams.update({'font.size': 22})
distance_to_solidangle = {10: 2.5183, 30: 1.4871, 50: 0.8111, 100: 0.3238, 200: 0.0881}
solidangle_values = list(distance_to_solidangle.values())
distance_values = list(distance_to_solidangle.keys())

# Plot FEPEint vs solid angle (Sigma) above
plt.figure()
plt.subplot(2, 1, 1)
for i, dist in enumerate(np.unique(distance)):
    mask = distance == dist
    x = solidangle[mask]
    y = fepeint[mask]
    plt.plot(x, y, marker='o', color='blue', label=f'{dist} mm')
    # plot linear fit
    popt, _ = curve_fit(linear_fit, x, y)
    linspace = np.linspace(min(x), max(x), 100)
    plt.plot(linspace, linear_fit(linspace, *popt), color=colors[i], linestyle='--', linewidth=1)

plt.xlabel('\u03A9 (sr)')
plt.ylabel(r'$\mathbf{\epsilon}_{\mathbf{FEP,int}}$')
plt.tick_params(axis='x', bottom=False, top=True, labelbottom=False, labeltop=True)
plt.xticks(solidangle_values, solidangle_values)
plt.grid(True)

# Move the x-axis label above the top x ticks
plt.gca().xaxis.set_label_coords(0.5, 1.16)

# Plot FEPE vs distance below
plt.subplot(2, 1, 2)
for i, dist in enumerate(np.unique(distance)):
    mask = distance == dist
    x = distance[mask]
    y = fepe[mask]
    plt.plot(x, y, marker='o', color='blue', label=f'{dist} mm')
    # plot log fit
    popt, _ = curve_fit(log_fit, x, y)
    linspace = np.linspace(min(x), max(x), 100)
    plt.plot(linspace, log_fit(linspace, *popt), color=colors[i], linestyle='--', linewidth=1)

plt.xlabel('Distance (mm)')
plt.ylabel(r'$\mathbf{\epsilon}_{\mathbf{FEP}}$')
plt.show()

""" # Plot FEPEint vs solid angle (Sigma) above
for i, en in enumerate(np.unique(energy)):
    plt.figure()
    plt.subplot(2, 1, 1)
    mask = energy == en
    x = solidangle[mask]
    y = fepeint[mask]
    plt.plot(x, y, marker = 'o', color= 'blue', label=f'{en} keV')
    linspace = np.linspace(min(x), max(x), 100)
    popt, _ = curve_fit(linear_fit, x, y)
    plt.plot(linspace, linear_fit(linspace, *popt), color='blue', linestyle='--', linewidth=1)

    plt.xlabel('\u03A9 (sr)')
    plt.ylabel(r'$\mathbf{\epsilon}_{\mathbf{FEP,int}}$')
    plt.tick_params(axis='x', bottom=False, top=True, labelbottom=False, labeltop=True)
    plt.xticks(solidangle_values, solidangle_values)
    plt.grid(True)

    plt.gca().xaxis.set_label_coords(0.5, 1.16)

    # Plot FEPE vs distance below
    plt.subplot(2, 1, 2)
    mask = energy == en
    x = distance[mask]
    y = fepe[mask]
    plt.plot(x, y, marker = 'o', color= 'blue', label=f'{en} keV')
    linspace = np.linspace(min(x), max(x), 100)
    popt, _ = curve_fit(log_fit, x, y)
    plt.plot(linspace, log_fit(linspace, *popt), color='blue', linestyle='--', linewidth=1)

    plt.xlabel('Distance (mm)')
    plt.ylabel(r'$\mathbf{\epsilon}_{\mathbf{FEP}}$')
    plt.grid(True)

    plt.tight_layout()
    plt.show() """



# Filter data for energies 443.98 and 1408
mask_443 = energy == 121.78
mask_1408 = energy == 1408

# Plot FEPEint vs solid angle (Sigma) above
plt.figure()
plt.subplot(2, 1, 1)
x_443 = solidangle[mask_443]
y_443 = fepeint[mask_443]
plt.plot(x_443, y_443, markers[0], color=colours[0], label='121.78 keV')
popt, _ = curve_fit(linear_fit, x_443, y_443)
linspace = np.linspace(min(x_443), max(x_443), 100)
plt.plot(linspace, linear_fit(linspace, *popt), color=colours[0], linestyle='--', linewidth=1)

# Plot data for energy 1408
x_1408 = solidangle[mask_1408]
y_1408 = fepeint[mask_1408]
plt.plot(x_1408, y_1408, markers[1], color=colours[1], label='1408 keV')
popt, _ = curve_fit(linear_fit, x_1408, y_1408)
plt.plot(linspace, linear_fit(linspace, *popt), color=colours[1], linestyle='--', linewidth=1)

plt.xlabel('\u03A9 (sr)')
plt.ylabel(r'$\mathbf{\epsilon}_{\mathbf{FEP,int}}$')
plt.tick_params(axis='x', bottom=False, top=True, labelbottom=False, labeltop=True)
plt.xticks(solidangle_values, solidangle_values)
plt.grid(True)

# Move the x-axis label above the top x ticks
plt.gca().xaxis.set_label_coords(0.5, 1.2)

# Plot FEPE vs distance below
plt.subplot(2, 1, 2)

# Plot data for energy 443.98
x_443 = distance[mask_443]
y_443 = fepe[mask_443]
plt.plot(x_443, y_443, markers[0], color=colours[0], label='121.78 keV')
linspace = np.linspace(min(x_443), max(x_443), 100)
popt, _ = curve_fit(log_fit, x_443, y_443)
plt.plot(linspace, log_fit(linspace, *popt), color=colours[0], linestyle='--', linewidth=1)

# Plot data for energy 1408
x_1408 = distance[mask_1408]
y_1408 = fepe[mask_1408]
plt.plot(x_1408, y_1408, markers[1], color=colours[1], label='1408 keV')
popt, _ = curve_fit(log_fit, x_1408, y_1408)
plt.plot(linspace, log_fit(linspace, *popt), color=colours[1], linestyle='--', linewidth=1)

plt.xlabel('Distance (mm)')
plt.ylabel(r'$\mathbf{\epsilon}_{\mathbf{FEP}}$')
plt.grid(True)

plt.tight_layout()
#plot the legend outside the plot at the bottom
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2)
plt.show()
