# ----------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------

import pypar
import numpy	#For numpy array definition

def BroadcastOneDouble(double):
  buff = numpy.array(range(1)).astype('f')
  buff[0] = double
  pypar.broadcast(buff, 0)
  double = float(buff[0])

  return double
