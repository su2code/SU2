import numpy as np

def lhc_unif(XB,NS,XI=None,maxits=10):
    ''' XS = lhc_unif(XB,NS,XI=None,maxits=10):
        
        Latin Hypercube Sampling with uniform density
        Iterates to maximize minimum L2 distance
        Accepts an array of points to respect while sampling
        
        Inputs:
            XB          - ndim x 2 array of [lower,upper] bounds
            NS          - number of new points to sample
            XI = None   - ni x ndim array of initial points to respect
            maxits = 10 - maximum number of iterations
        
        Outputs:
            XS - ns x ndim array of sampled points
    '''
    
    # dimension
    XB = np.atleast_2d(XB)
    ND = XB.shape[0]
    
    # initial points to respect
    if XI is None:
        XI = np.empty([0,ND])
    else:
        XI = np.atleast_2d(XI)
       
    # output points
    XO = []
    
    # initialize
    mindiff = 0;
    
    # maximize minimum distance
    for it in range(maxits):
        
        # samples
        S = np.zeros([NS,ND])
        
        # populate samples
        for i_d in range(ND):
            S[:,i_d] = ( np.random.random([1,NS]) + np.random.permutation(NS) ) / NS
        XS = S*(XB[:,1]-XB[:,0]) + XB[:,0]        
        
        # add initial points
        XX = np.vstack([ XI , XS ])
        
        # calc distances
        vecdiff = vec_dist(XX)[0]
        
        # update
        if vecdiff > mindiff:
            mindiff = vecdiff
            XO = XX
        
    #: for iterate
    
    return XO

def vec_dist(X,P=None):
    ''' calculates distance between points in matrix X 
        with each other, or optionally to given point P
        returns min, max and matrix/vector of distances
    '''
    
    # distance matrix among X
    if P is None:
        
        nK,nD = X.shape
        
        d = np.zeros([nK,nK,nD])
        for iD in range(nD):
            d[:,:,iD] = np.array([X[:,iD]])-np.array([X[:,iD]]).T
        D = np.sqrt( np.sum( d**2 , 2 ) )
        
        diag_inf = np.diag( np.ones([nK])*np.inf )
        dmin = np.min(np.min( D + diag_inf ))
        dmax = np.max(np.max( D ))
        
    # distance vector to P
    else:
        assert P.shape[0] == 1 , 'P must be a horizontal vector'
        D = np.array([ np.sqrt( np.sum( (X-P)**2 , 1 ) ) ]).T
        dmin = D.min()
        dmax = D.max()
        
    return (dmin,dmax,D)