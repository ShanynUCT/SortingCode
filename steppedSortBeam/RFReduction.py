import numpy as np
import pandas as pd
from pathlib import Path
import uproot
import sys
import ROOT
import csv
from matplotlib import pyplot as plt

# Read in data
data = sys.argv[1]
save_dir = Path(data).parent / 'plots/RFReduction'

if save_dir.exists():
    print(f"\nThe directory {save_dir} already exists.")
else:
    save_dir.mkdir(parents=True, exist_ok=False)
    print(f"\nThe directory {save_dir} was created.")

f = uproot.open(data)
directory = f["dir_aligned"]
tree = directory["AlignedData"]

alignedData = tree.arrays(['detectorID', 'globalTime', 'timeF', 'timeS', 'energyS'], library="pd")
alignedData = alignedData.drop_duplicates()
alignedData = alignedData.reset_index(drop=True)
#keep only detector 3 and 4
alignedData = alignedData[(alignedData['detectorID'] == 3) | (alignedData['detectorID'] == 4)]
alignedData['RFtime'] = np.nan
alignedData['RFtime'] = alignedData.loc[alignedData['detectorID'] == 4, 'globalTime']
alignedData['LaBr3time'] = alignedData.loc[alignedData['detectorID'] == 3, 'globalTime']
alignedData['RFtime'] = alignedData['RFtime'].bfill()
# remove duplicate rows

# drop timeF, timeS, and globalTime
alignedData = alignedData.drop(columns=['timeF', 'timeS', 'globalTime'])
alignedData = alignedData[alignedData['detectorID'] != 4]
alignedData = alignedData[(alignedData['energyS'] != 0) | (alignedData['energyS'].isna() == False)]

alignedData['timediff43'] = alignedData['RFtime'] - alignedData['LaBr3time']

print(alignedData)

# Plotting
plt.figure()
plt.hist(alignedData['timediff43'], bins=30000, range=(0, 300000))
plt.xlabel('Time Difference (ns)')
plt.ylabel('Counts')
plt.title('Time Difference Between RF and LaBr3')
plt.yscale('log')
plt.show()