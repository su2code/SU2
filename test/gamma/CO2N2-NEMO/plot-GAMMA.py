import sys
import numpy as np
import itertools as IT
import matplotlib
from matplotlib import rc, font_manager
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams.update({'font.size': 19})
import matplotlib.ticker
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import pandas as pd
import pathlib
import matplotlib.patches as patches

fname = 'gamma.csv'
dataset = pd.read_csv(fname)
params_value = pd.DataFrame(dataset)
param1 = params_value.iloc[:,0]
param2 = params_value.iloc[:,1]

lines = plt.plot(param1, param2)
plt.setp(lines, color='b', linewidth=1.0, linestyle='-')

#name = 'gamma-ref.csv'
#ataset = pd.read_csv(fname)
#arams_value = pd.DataFrame(dataset)
#aram1 = params_value.iloc[:,0]
#aram2 = params_value.iloc[:,1]

#ines = plt.plot(param1, param2, marker="o")

#fname = '../gamma.csv'
#dataset = pd.read_csv(fname)
#params_value = pd.DataFrame(dataset)
#param1 = params_value.iloc[:,0]
#param2 = params_value.iloc[:,1]
#
#lines = plt.plot(param1, param2)
#plt.setp(lines, color='r', linewidth=1.0, linestyle='-')






area = 50
rc('text.latex', preamble=r'\usepackage{cmbright}')
plt.rc('grid', linestyle=":", color='grey')
plt.grid(False)

plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)

t1 = mlines.Line2D([], [], color='b', linestyle='-', label='SU2-MPP')
t2 = mlines.Line2D([], [], color='r', linestyle='-', label='MPP')
#t3 = mlines.Line2D([], [], color='orange', marker="o", label='ref table')
#t4 = mlines.Line2D([], [], color='k', linestyle='-', label='Level 3')
plt.legend(handles=[t1,t2], frameon=False) 

plt.ylabel(r'$\gamma$', fontsize=30)
plt.xlabel(r'$T$ [K]', fontsize=30)
plt.xlim(0, 20000)

plt.show()
print('')