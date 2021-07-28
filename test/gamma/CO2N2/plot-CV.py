import sys
import numpy as np
import itertools as IT
import matplotlib
from matplotlib import rc, font_manager
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams.update({'font.size': 30})
import matplotlib.ticker
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import pandas as pd
import pathlib
import matplotlib.patches as patches

fname = 'cv-su2-pig.csv'
dataset = pd.read_csv(fname)
params_value = pd.DataFrame(dataset)
param1 = params_value.iloc[:,0]
param2 = params_value.iloc[:,1]

lines = plt.plot(param1, param2)
plt.setp(lines, color='gray', linewidth=2.0, linestyle='-')

fname = 'cv-su2.csv'
dataset = pd.read_csv(fname)
params_value = pd.DataFrame(dataset)
param1 = params_value.iloc[:,0]
param2 = params_value.iloc[:,1]

lines = plt.plot(param1, param2)
plt.setp(lines, color='k', linewidth=2.0, linestyle='-')



fname = 'cv-su2-nemo.csv'
dataset = pd.read_csv(fname)
params_value = pd.DataFrame(dataset)
param1 = params_value.iloc[:,0]
param2 = params_value.iloc[:,1]
param3 = params_value.iloc[:,2]

lines = plt.plot(param1, param2)
plt.setp(lines, color='r', linewidth=2.0, linestyle='-')

lines = plt.plot(param1, param3)
plt.setp(lines, color='b', linewidth=2.0, linestyle='-')


#param4 = params_value.iloc[:,1] + params_value.iloc[:,2]
#lines = plt.plot(param1, param4)
#plt.setp(lines, color='k', linewidth=2.0, linestyle='-', marker= 'o' )





area = 50
rc('text.latex', preamble=r'\usepackage{cmbright}')
plt.rc('grid', linestyle=":", color='grey')
plt.grid(True)

plt.gca().spines['right'].set_visible(True)
plt.gca().spines['top'].set_visible(True)

t1 = mlines.Line2D([], [], color='gray', linestyle='-', label='Cv PIG')
t2 = mlines.Line2D([], [], color='k', linestyle='-', label='Cv TPG')
t3 = mlines.Line2D([], [], color='r', linestyle='-', label='Cv$_{tr}$ NEQ')
t4 = mlines.Line2D([], [], color='b', linestyle='-', label='Cv$_{ve}$ NEQ')

#t4 = mlines.Line2D([], [], color='k', linestyle='-', label='Level 3')
plt.legend(handles=[t1,t2,t3,t4], frameon=False) 

plt.ylabel(r'C$_v$ [J/(kg K)]' , fontsize=30)
plt.xlabel(r'$T$ [K]', fontsize=30)
plt.xlim(100, 10000)
plt.ylim(0, 1800)

plt.show()
print('')