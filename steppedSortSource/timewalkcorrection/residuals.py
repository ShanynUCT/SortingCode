import matplotlib.pyplot as plt
import numpy as np


plt.rcParams["axes.labelweight"] = "bold"
# Make all plots square
#make all fonts 22
plt.rcParams.update({'font.size': 28})
plt.tick_params(axis='both', which='major', labelsize=26)
plt.tick_params(axis='both', which='major', labelsize=26)

plt.rcParams['figure.figsize'] = [10, 10]

# Data for ch1slowL2
energy_ch1slowL2 = np.array([41.1, 121.8, 244.7, 295.9, 344.3, 411, 444, 778.9, 867.4, 964.1, 1099, 1408])
pol1_ch1slowL2 = np.array([-0.45, 1.23, -0.29, 1.25, -1.71, -2.57, -2.06, -2.70, 0.36, -1.49, -3.03, -0.20])
pol2_ch1slowL2 = np.array([-0.22, 1.86, 0.82, 2.52, -0.32, -1.06, -0.52, -1.39, 1.41, -0.81, -3.07, -2.55])
pol3_ch1slowL2 = np.array([-0.33, 1.62, 0.49, 2.19, -0.63, -1.31, -0.74, -1.15, 1.75, -0.42, -2.72, -3.05])

# Data for ch0slowL0
energy_ch0slowL0 = np.array([41.1, 121.8, 244.7, 295.9, 344.3, 411, 444, 778.9, 867.4, 964.1, 1099, 1408])
pol1_ch0slowL0 = np.array([-0.32, 0.06, -1.21, -2.11, -2.18, -4.17, -2.69, -2.05, -1.08, -1.61, -2.35, 1.04])
pol2_ch0slowL0 = np.array([0.08, 1.14, 0.70, 0.07, 0.19, -1.59, -0.04, 0.19, 0.71, -0.48, -2.42, -2.98])
pol3_ch0slowL0 = np.array([0.11, 1.22, 0.81, 0.17, 0.29, -1.50, 0.03, 0.11, 0.60, -0.61, -2.54, -2.82])

# Data for ch8fastL0
energy_ch8fastL0 = np.array([41.1, 121.8, 244.7, 344.3, 411, 444, 778.9, 867.4, 964.1, 1099, 1408])
pol1_ch8fastL0 = np.array([-3.12, -1.01, -1.58, -2.91, -6.74, -2.75, -1.76, 0.17, -1.39, -3.47, -0.63])
pol2_ch8fastL0 = np.array([-2.64, 0.23, 0.60, -0.20, -3.79, 0.28, 0.84, 2.28, -0.02, -3.46, -5.07])
pol3_ch8fastL0 = np.array([-2.30, 1.00, 1.65, 0.79, -2.95, 1.03, 0.25, 1.41, -1.08, -4.45, -3.71])

# Data for ch9fastL2
energy_ch9fastL2 = np.array([41.1, 121.8, 244.7, 295.9, 344.3, 411, 444, 778.9, 867.4, 964.1, 1099, 1408])
pol1_ch9fastL2 = np.array([-1.00, 0.57, -0.44, 1.55, -1.65, -2.98, -1.24, -3.03, -0.58, -1.60, -2.26, 0.99])
pol2_ch9fastL2 = np.array([-0.73, 1.28, 0.81, 2.96, -0.10, -1.29, 0.49, -1.56, 0.60, -0.86, -2.30, -1.63])
pol3_ch9fastL2 = np.array([-0.90, 0.89, 0.28, 2.44, -0.58, -1.69, 0.14, -1.18, 1.12, -0.24, -1.74, -2.41])

# Create a list of data sets for convenience
datasets = [
    ("ch1slowL2", energy_ch1slowL2, pol1_ch1slowL2, pol2_ch1slowL2, pol3_ch1slowL2),
    ("ch0slowL0", energy_ch0slowL0, pol1_ch0slowL0, pol2_ch0slowL0, pol3_ch0slowL0),
    ("ch8fastL0", energy_ch8fastL0, pol1_ch8fastL0, pol2_ch8fastL0, pol3_ch8fastL0),
    ("ch9fastL2", energy_ch9fastL2, pol1_ch9fastL2, pol2_ch9fastL2, pol3_ch9fastL2)
]

# Create separate plots for each dataset
for name, energy, pol1, pol2, pol3 in datasets:
    plt.figure(figsize=(10, 6))
    plt.scatter(energy, pol1, color='blue', marker='o', label='y = ax + b')
    plt.scatter(energy, pol2, color='green', marker='s', label='y = ax$^2$ + bx + c')
    plt.scatter(energy, pol3, color='red', marker='^', label='y = ax$^3$ + bx$^2$ + cx + d')
    plt.axhline(0, color='black', linewidth=2)  # Dark black line at y=0
    plt.xlabel('Energy (keV)')
    plt.ylabel('Residuals (keV)')
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=3, fontsize=16)
    plt.grid(True)
    plt.tight_layout()
    plt.show()
