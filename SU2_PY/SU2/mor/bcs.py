import random as rand
from math import floor, log10

class isothermal_bc:
    """Usage: 
        samples = isothermal_bc(mu_high, mu_low, nparams, nsamples)
    Output:
        list of nsamples (random) sample strings, each string of length nparams
    """
    def __init__(self, mu_high, mu_low, nparams, nsamples):
        self.mu_high = mu_high
        self.mu_low = mu_low
        self.nparams = nparams
        self.nsamples = nsamples
        
        
        
    def get_samples_list(self):
        """Returns list of samples
        """
        samples = []
        for i in range(self.nsamples):
            samples.append(self.__get_sample())
        
        return samples
        
    def __get_sample(self):
        """Creates one sample of length nparams with random values
           Assumes each parameter has the same max and min values
        """
        #d = "MARKER_ISOTHERMAL= ( " #wall1, 5.0, wall2, 10.0, wall3, 5.0, wall4, 10.0 ) "
        d = "( "
        def round_sig(x, sig=2):
            return round(x, sig-int(floor(log10(abs(x))))-1)
            
        for i in range(self.nparams):
            val = rand.random() * (self.mu_high-self.mu_low) + self.mu_low
            d += "wall" + str(i+1) + ", " + str(round_sig(val)) + ", "
            
        d = d[:-2] + " )"
        
        return d
        