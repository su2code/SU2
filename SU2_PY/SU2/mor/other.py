class DoE:
    def __init__(self):
        pass
        
    def lhs_samples(self, x=3, n=3, xlim=[0.0,0.1]):
        """
        INPUTS:
            x: dimensions (number of parameters)
            n: factoring for each dimension 
        """
        #for i in range(x):
        pass
        
        
class shape_params: 
    """
    Helper class for creation of hicks-henne bump functions
    """
    def __init__(self, up_loc=0.5, low_loc=0.5):
        
        self.DV_PARAM = self.__dv_params(up_loc, low_loc)
        self.n = len(up_loc) + len(low_loc)
        self.DV_KIND = self.__dv_kind()
        
    def lhs_samples(self, x=3, n=3, xlim=[0.0,0.1]):
        """
        INPUTS:
            x: dimensions (number of parameters)
            n: factoring for each dimension 
        """
        pass

    def __dv_params(self, up_loc, low_loc):
        '''*_vals: array of hicks henne locations on * surface'''
        
        # see SU2/io/config.py line 183
        d = {}
        d['PARAM'] = []
        d['FFDTAG'] = []
        d['SIZE'] = []
        for loc in up_loc:
            d['PARAM'].append([1.0,loc])
            d['FFDTAG'].append([])
            d['SIZE'].append(1)
            
        for loc in low_loc:
            d['PARAM'].append([0.0,loc])
            d['FFDTAG'].append([])
            d['SIZE'].append(1)
            
        return d
        
    def get_dv_values(self, up_limit=0.05):
        '''up_limit: upper limit of possible hicks henne bump magnitude
           lower limit assumed to be zero
           magnitudes are found at random
        '''
        val = []
        for i in range(self.n):
            val[i] = rand.random()*up_limit
        
        d = ''
        for v in val:
            d += str(val) + ', ' 
            
        return d[0:-2]
        
    def __dv_kind(self):
        '''n: number of hicks henne bumps to use'''
        h = 'HICKS_HENNE'
        d = ''
        for i in range(self.n):
            d += h + ', '
            
        return d[0:-2]

    def get_data(self):
        return self.DV_PARAM, self.DV_KIND