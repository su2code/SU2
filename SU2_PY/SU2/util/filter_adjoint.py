#!/usr/bin/env python

## \file filter_adjoint.py
#  \brief Applies various filters to the adjoint surface sensitivities of an airfoil
#  \author T. Lukaczyk, F. Palacios
#  \version 7.5.1 "Blackbird"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2023, SU2 Contributors (cf. AUTHORS.md)
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.

from __future__ import division, print_function, absolute_import

import os
import math
from numpy import pi
from optparse import OptionParser
import numpy as np
import libSU2
import libSU2_mesh

# plotting with matplotlib
try:
    import pylab as plt
    pylab_imported = True
except ImportError:
    pylab_imported = False


# -------------------------------------------------------------------
#  MAIN
# -------------------------------------------------------------------
def main():

    # Command Line Options
    parser=OptionParser()
    parser.add_option( "-f", "--file",       dest="filename",
                       help="read config from FILE", metavar="FILE" )
    parser.add_option( "-t", "--type",       dest="filter_type",  default='LAPLACE',
                       help="apply filter TYPE", metavar="TYPE" )
    parser.add_option( "-m", "--marker",     dest="marker_name",  default='airfoil',
                        help="use marker named TAG", metavar="TAG" )
    parser.add_option( "-c", "--chord",      dest="chord_length", default=1.0,
                        help="reference CHORD length", metavar="CHORD" )

    (options, args)=parser.parse_args()
    options.chord_length = float( options.chord_length )

    # run filter
    process_surface_adjoint( options.filename     ,
                             options.filter_type  ,
                             options.marker_name  ,
                             options.chord_length  )

#: def main()


# -------------------------------------------------------------------
#  PROCESS SURFACE ADJOINT
# -------------------------------------------------------------------
def process_surface_adjoint( config_filename       ,
                             filter_type='LAPLACE' ,
                             marker_name='airfoil' ,
                             chord_length=1.0       ):

    print('')
    print('-------------------------------------------------------------------------')
    print('|                 SU2 Suite (Process Surface Adjoint)                   |')
    print('-------------------------------------------------------------------------')
    print('')

    # some other defaults
    c_clip   = 0.01 # percent chord to truncate
    fft_copy = 5    # number of times to copy the fft signal
    smth_len = 0.05 # percent chord smoothing window length
    lapl_len = 1e-4 # laplace smoothing parameter

    # read config file
    config_data = libSU2.Get_ConfigParams(config_filename)
    surface_filename = config_data['SURFACE_ADJ_FILENAME'] + '.csv'
    print(surface_filename)
    mesh_filename    = config_data['MESH_FILENAME']
    gradient         = config_data['OBJECTIVE_FUNCTION']

    print('Config filename = %s' % config_filename)
    print('Surface filename = %s' % surface_filename)
    print('Filter Type = %s' % filter_type)

    # read adjoint data
    adj_data = np.genfromtxt( surface_filename    ,
                              dtype       = float ,
                              delimiter   = ','   ,
                              skip_header = 1      )

    # read mesh data
    mesh_data = libSU2_mesh.Read_Mesh(mesh_filename)

    # proces adjoint data
    P      = map(int, adj_data[:,0] )
    X      = adj_data[:,6].copy()
    Y      = adj_data[:,7].copy()
    Sens   = adj_data[:,1].copy()
    PsiRho = adj_data[:,2].copy()
    I      = range(0,len(P)) # important - for unsorting durring write

    # store in dict by point index
    adj_data_dict = dict( zip( P , zip(X,Y,Sens,PsiRho,I) ) )

    # sort airfoil points
    iP_sorted,_ = libSU2_mesh.sort_Airfoil(mesh_data,marker_name)
    assert(len(iP_sorted) == len(P))

    # rebuild airfoil loop
    i = 0
    for this_P in iP_sorted:
        # the adjoint data entry
        this_adj_data = adj_data_dict[this_P]
        # re-sort
        P[i]      = this_P
        X[i]      = this_adj_data[0]
        Y[i]      = this_adj_data[1]
        Sens[i]   = this_adj_data[2]
        PsiRho[i] = this_adj_data[3]
        I[i]      = this_adj_data[4]
        # next
        i = i+1
    #: for each point

    # calculate arc length
    S = np.sqrt( np.diff(X)**2 + np.diff(Y)**2 ) / chord_length
    S = np.cumsum( np.hstack([ 0 , S ]) )

    # tail trucating, by arc length
    I_clip_lo = S < S[0]  + c_clip
    I_clip_hi = S > S[-1] - c_clip
    S_clip    = S.copy()
    Sens_clip = Sens.copy()
    Sens_clip[I_clip_hi] = Sens_clip[I_clip_hi][0]
    Sens_clip[I_clip_lo] = Sens_clip[I_clip_lo][-1]


    # some edge length statistics
    dS_clip = np.diff(S_clip)
    min_dS  = np.min ( dS_clip )
    mean_dS = np.mean( dS_clip )
    max_dS  = np.max ( dS_clip )
    #print 'min_dS = %.4e ; mean_dS = %.4e ; max_dS = %.4e' % ( min_dS , mean_dS , max_dS )

    # --------------------------------------------
    #  APPLY FILTER

    if filter_type == 'FOURIER':
        Freq_notch = [ 1/max_dS, np.inf ] # the notch frequencies
        Sens_filter,Frequency,Power = fft_filter( S_clip,Sens_clip, Freq_notch, fft_copy )
        #Sens_filter = smooth(S_clip,Sens_filter, 0.03,'blackman') # post smoothing

    elif filter_type == 'WINDOW':
        Sens_filter   = window( S_clip, Sens_clip, smth_len, 'blackman' )

    elif filter_type == 'LAPLACE':
        Sens_filter   = laplace( S_clip, Sens_clip, lapl_len )

    elif filter_type == 'SHARPEN':
        Sens_smooth   = smooth( S_clip, Sens_clip  , smth_len/5, 'blackman' ) # pre smoothing
        Sens_smoother = smooth( S_clip, Sens_smooth, smth_len  , 'blackman' )
        Sens_filter = Sens_smooth + (Sens_smooth - Sens_smoother)             # sharpener
    else:
        raise Exception('unknown filter type')

    # --------------------------------------------
    #  PLOTTING

    if pylab_imported:

        # start plot
        fig = plt.figure(gradient)
        plt.clf()
        #if not fig.axes:          # for comparing two filter calls
            #plt.subplot(1,1,1)
        #ax = fig.axes[0]
        #if len(ax.lines) == 4:
            #ax.lines.pop(0)
            #ax.lines.pop(0)

        # SENSITIVITY
        plt.plot(S     ,Sens       ,color='b') # original
        plt.plot(S_clip,Sens_filter,color='r') # filtered

        plt.xlim(-0.1,2.1)
        plt.ylim(-5,5)
        plt.xlabel('Arc Length')
        plt.ylabel('Surface Sensitivity')

        #if len(ax.lines) == 4:
            #seq = [2, 2, 7, 2]
            #ax.lines[0].set_dashes(seq)
            #ax.lines[1].set_dashes(seq)

        plot_filename = os.path.splitext(surface_filename)[0] + '.png'
        plt.savefig('Sens_'+plot_filename,dpi=300)

        # zoom in
        plt.ylim(-0.4,0.4)
        plt.savefig('Sens_zoom_'+plot_filename,dpi=300)

        # SPECTRAL
        if filter_type == 'FOURIER':

            plt.figure('SPECTRAL')
            plt.clf()

            plt.plot(Frequency,Power)

            #plt.xlim(0,Freq_notch[0]+10)
            plt.xlim(0,200)
            plt.ylim(0,0.15)

            plt.xlabel('Frequency (1/C)')
            plt.ylabel('Surface Sensitivity Spectal Power')

            plt.savefig('Spectral_'+plot_filename,dpi=300)

        #: if spectral plot

    #: if plot

    # --------------------------------------------
    #  SAVE SURFACE FILE

    # reorder back to input surface points
    Sens_out    = np.zeros(len(S))
    Sens_out[I] = Sens_filter # left over from sort
    adj_data[:,1] = Sens_out

    # get surface header
    surface_orig = open(surface_filename,'r')
    header = surface_orig.readline()
    surface_orig.close()

    # get list of prefix names
    prefix_names = libSU2.get_AdjointPrefix(None)
    prefix_names = prefix_names.values()

    # add filter prefix, before adjoint prefix
    surface_filename_split = surface_filename.rstrip('.csv').split('_')
    if surface_filename_split[-1] in prefix_names:
        surface_filename_split = surface_filename_split[0:-1] + ['filtered'] + [surface_filename_split[-1]]
    else:
        surface_filename_split = surface_filename_split + ['filtered']
    surface_filename_new = '_'.join(surface_filename_split) + '.csv'

    # write filtered surface file (only updates Sensitivity)
    surface_new = open(surface_filename_new,'w')
    surface_new.write(header)
    for row in adj_data:
        for i,value in enumerate(row):
            if i > 0:
                surface_new.write(', ')
            if i == 0:
                surface_new.write('%i' % value )
            else:
                surface_new.write('%.16e' % value )
        surface_new.write('\n')
    surface_new.close()


    print('')
    print('----------------- Exit Success (Process Surface Adjoint) ----------------')
    print('')

#: def process_surface_adjoint()


# -------------------------------------------------------------------
#  LAPLACIAN SMOOTHING
# -------------------------------------------------------------------

def laplace(t,x,e):
    ''' Laplacian filter
        input:
            t - time sample vector
            x - signal vector x(t)
            e - smoother coefficient (e>0)

        output:
            y: smoothed signal at t
    '''

    n_x = len(x)

    # padding
    t_1 = t[ 0] + t[-2]-t[-1]
    t_2 = t[-1] + t[ 1]-t[ 0]
    t_p = np.hstack([ t_1  , t , t_2   ])
    x_p = np.hstack([ x[0] , x , x[-1] ])

    # finite differencing
    dt_f = t_p[2:  ] - t_p[1:-1]
    dt_b = t_p[1:-1] - t_p[0:-2]
    dt_c = t_p[2:  ] - t_p[0:-2]

    # diagonal coefficients
    Coeff = e * 2.0 / (dt_b*dt_f*dt_c)
    diag_c =  Coeff*dt_c
    diag_f = -Coeff*dt_b
    diag_b = -Coeff*dt_f

    # system matrix
    A = ( np.diag(diag_c      , 0) +
          np.diag(diag_f[0:-1], 1) +
          np.diag(diag_b[1:  ],-1) +
          np.diag(np.ones(n_x), 0)  )

    # periodic conditions
    #A[1,-1] = dt_b[0]
    #A[-1,1] = dt_f[-1]

    # rhs
    b = np.array([x]).T

    # boundary conditions

    # signal start
    i_d = 0
    A[i_d,:] = 0.0
    A[i_d,i_d]   = 1.0        # dirichlet
    #A[i_d,i_d+1] =  1.0       # neuman
    #A[i_d,i_d]   = -1.0
    #b[i_d] = 0.0 #x[i_d+1]-x[i_d]

    # signal end
    i_d = n_x-1
    A[i_d,:] = 0.0
    A[i_d,i_d]   = 1.0        # dirichlet
    #A[i_d,i_d]   =  1.0       # neuman
    #A[i_d,i_d-1] = -1.0
    #b[i_d] = 0.0 #x[i_d]-x[i_d-1]

    # solve
    y = np.linalg.solve(A,b)
    y = y[:,0]

    return y

#: def laplace


# -------------------------------------------------------------------
#  FFT NOTCH FILTER
# -------------------------------------------------------------------

def fft_filter(t,x,n,c=1):
    ''' Notch filter with Fast Fourier Transform
        input:
            t = input time vector
            x = input signal vector
            n = [low,high] frequency range to supress
            c = number of times to duplicate signal

        output:
            y = smoothed signal at t

        signal will be interpolated to constant spacing
    '''

    assert(c>0)

    # choose sampling frequency
    min_dt  = np.min( np.diff(t) )
    Ts = min_dt/2
    Fs = 1/Ts

    # interpolate to constant spacing
    nt_lin = int( t[-1]/Ts )
    t_lin = np.linspace(0,t[-1],nt_lin)
    x_lin = np.interp(t_lin,t,x)

    # pad last index
    t_lin = np.hstack([ t_lin , t_lin[0:10]+t_lin[-1] ])
    x_lin = np.hstack([ x_lin , np.ones(10)*x_lin[-1] ])

    # copy signal
    for ic in range(c-1):
        t_lin = np.hstack([ t_lin[0:-2] , t_lin[1:]+t_lin[-1] ])
        x_lin = np.hstack([ x_lin[0:-2] , x_lin[1:] ])
    nt = len(t_lin)

    # perform fourier transform
    nxtpow2 = int(math.log(nt, 2))+1 # next power of 2
    nfft = 2**nxtpow2                # fft efficiency
    P = np.fft.rfft(x_lin,nfft)      # the transform
    a = np.angle(P)                  # complex
    p = np.absolute(P)               # complex
    p = p/nt                         # normalize
    p = 2*p[0:(nfft/2)]              # symmetric
    a = a[0:(nfft/2)]                # symmetric

    # frequency domain
    F = np.arange(0,nfft/2) * Fs/nfft

    # for return
    Freq = F.copy()

    # --------------------------------------------------
    #  THE NOTCH FILTER

    # filter multiplier
    k = np.ones(nfft/2)

    # clip power within notch frequencies
    I_fil = np.logical_and( F>n[0] , F<n[1] )
    k[I_fil] = 0.0

    # change the power spectrum
    p = p*k

    # For Return
    Pow = p.copy()

    # untransform
    p = p*nt/2.
    p = np.hstack( [ p , p[::-1] ] )
    a = np.hstack( [ a , -a[::-1] ] )
    P = p*(np.cos(a) + 1j*np.sin(a))
    y_lin = np.fft.irfft(P,nfft)       # the inverse transform
    y_lin = y_lin[0:nt]

    # interpolate back to given t
    y = np.interp(t,t_lin,y_lin)

    return y,Freq,Pow

# def: fft_filter()


# -------------------------------------------------------------------
#  WINDOWED SMOOTHING
# -------------------------------------------------------------------

def window(t,x,window_delta,window='hanning'):
    """Smooth the data using a window with requested size and shape

    original source:
    http://www.scipy.org/Cookbook/SignalSmooth

    input:
        t: input time samples
        x: input signal at t
        window_delta: length (in units of t) of the window
        window: type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        y: smoothed signal at t

    """

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if window not in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is not of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

    # interpolate to constant time sample width
    min_dt = np.min( np.diff(t) )
    Ts = min_dt/2
    nt_lin = int( t[-1]/Ts )
    t_lin = np.linspace(0,t[-1],nt_lin)
    x_lin = np.interp(t_lin,t,x)

    # window sample length
    window_len = int( window_delta / Ts )

    # padding
    s=np.r_[x_lin[window_len-1:0:-1],x_lin,x_lin[-1:-window_len:-1]]

    # window template
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    # the filter
    y_lin = np.convolve(w/w.sum(),s,mode='valid')

    # remove padding
    y_lin = y_lin[((window_len-1)/2):-(window_len/2)]

    # interpolate back to given t
    y = np.interp(t,t_lin,y_lin)

    return y

#: def window()


# -----------------------------------------------------------------
#  Run Main from Command Line
# -----------------------------------------------------------------

if __name__ == '__main__':
    main()


