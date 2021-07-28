import sys
import numpy as np
import itertools as IT
import matplotlib
from matplotlib import rc, font_manager
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams.update({'font.size': 35})
import matplotlib.ticker
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import pandas as pd
import pathlib
import matplotlib.patches as patches

fname = 'gamma-su2.csv'
dataset = pd.read_csv(fname)
params_value = pd.DataFrame(dataset)
param1 = params_value.iloc[:,0]
param2 = params_value.iloc[:,1]

fig,ax = plt.subplots()
ax.plot(param1, param2, color='k', linewidth=2.0, linestyle='-')

ax2=ax.twinx()

fname = 'cv-su2.csv'
dataset = pd.read_csv(fname)
params_value = pd.DataFrame(dataset)
param1 = params_value.iloc[:,0]
param2 = params_value.iloc[:,1]/190.258

ax2.plot(param1, param2, color='b', linewidth=2.0, linestyle='-')


area = 50
rc('text.latex', preamble=r'\usepackage{cmbright}')
plt.rc('grid', linestyle=":", color='grey')
plt.grid(False)

plt.gca().spines['right'].set_visible(True)
plt.gca().spines['top'].set_visible(True)

t1 = mlines.Line2D([], [], color='k', linestyle='-', label='$\gamma$')
t2 = mlines.Line2D([], [], color='b', linestyle='-', label='C$_v$')
##t3 = mlines.Line2D([], [], color='orange', marker="o", label='ref table')
#t4 = mlines.Line2D([], [], color='k', linestyle='-', label='Level 3')
plt.legend(handles=[t1,t2], frameon=False) 

ax.set_ylabel(r'$\gamma$', fontsize=35)
ax2.set_ylabel(r' C$_v$/R' , fontsize=35)
ax.set_xlabel(r'$T$ [K]', fontsize=35)
ax.set_xlim(100, 10000)

plt.show()
print('')