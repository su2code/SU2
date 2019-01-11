# test using navpack pod

from read_flow import Read_Flow
import numpy as np
import pylab
from navpack import dimred, snapshot


designs = 20
snapshots = list()
for m in range(designs):
    folder = 'DESIGNS/DSN_0'
    
    m += 1
    if m < 10:
        folder += '0' + str(m) + '/flow.dat'
    else:
        folder += str(m) + '/flow.dat'
        
    results = Read_Flow(folder,delimiter='","')
    
    # grab temperature data only, must be in numpy array format
    temp = np.asarray(results.return_data('Conservative_1'))
    snapshots.append(snapshot.Snapshot({"d": temp}))

# perform POD using Navpack pod routine    
pod = dimred.POD.CreateFromSnapshots(snapshots,"d")

x = np.asarray(results.return_data('x'))
y = np.asarray(results.return_data('y'))
modes = 8
for m in range(modes):
    pylab.subplot(2,4,m+1)
    pylab.scatter(x, y, c=pod.Phi[m])

pylab.suptitle("First 8 POD modes")
pylab.show()
    
    