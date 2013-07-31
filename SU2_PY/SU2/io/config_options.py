
from ..util import ordered_bunch

class OptionError(Exception):
    pass

class Option(object):
    
    def __init__(self):
        self.val = ""

    def __get__(self):
        return self.val

    def __set__(self,newval):
        self.val = newval

#: class Option

class MathProblem(Option):

    def __init__(self,*args,**kwarg):
        super(MathProblem,self).__init__(*args,**kwarg)
        self.validoptions = ['DIRECT','ADJOINT','LINEARIZED']

    def __set__(self,newval):
        if not self.newval in self.validoptions:
            raise OptionError("Invalid option. Valid options are: %s"%self.validoptions)
        super(MathProblem,self).__set__(newval)

#: class MathProblem



class DEFINITION_DV(ordered_bunch):
    """ SU2.io.config.DEFINITION_DV()
    
        List of design variables (Design variables are separated by semicolons)
        - HICKS_HENNE ( 1, Scale | Mark. List | Lower(0)/Upper(1) side, x_Loc )
        - COSINE_BUMP ( 2, Scale | Mark. List | Lower(0)/Upper(1) side, x_Loc, x_Size )
        - SPHERICAL ( 3, Scale | Mark. List | ControlPoint_Index, Theta_Disp, R_Disp )
        - NACA_4DIGITS ( 4, Scale | Mark. List |  1st digit, 2nd digit, 3rd and 4th digit )
        - DISPLACEMENT ( 5, Scale | Mark. List | x_Disp, y_Disp, z_Disp )
        - ROTATION ( 6, Scale | Mark. List | x_Axis, y_Axis, z_Axis, x_Turn, y_Turn, z_Turn )
        - FFD_CONTROL_POINT ( 7, Scale | Mark. List | Chunk, i_Ind, j_Ind, k_Ind, x_Mov, y_Mov, z_Mov )
        - FFD_DIHEDRAL_ANGLE ( 8, Scale | Mark. List | Chunk, x_Orig, y_Orig, z_Orig, x_End, y_End, z_End )
        - FFD_TWIST_ANGLE ( 9, Scale | Mark. List | Chunk, x_Orig, y_Orig, z_Orig, x_End, y_End, z_End )
        - FFD_ROTATION ( 10, Scale | Mark. List | Chunk, x_Orig, y_Orig, z_Orig, x_End, y_End, z_End )
        - FFD_CAMBER ( 11, Scale | Mark. List | Chunk, i_Ind, j_Ind )
        - FFD_THICKNESS ( 12, Scale | Mark. List | Chunk, i_Ind, j_Ind )
        - FFD_VOLUME ( 13, Scale | Mark. List | Chunk, i_Ind, j_Ind )
        - FOURIER ( 14, Scale | Mark. List | Lower(0)/Upper(1) side, index, cos(0)/sin(1) )
    """
    
    def __init__(self,*args,**kwarg):
        ordered_bunch.__init__(self)
        self.KIND   = []
        self.SCALE  = []
        self.MARKER = []
        self.PARAM  = []
        self.update(ordered_bunch(*args,**kwarg))
    
    def append(self,new_dv):
        self.KIND.  append(new_dv['KIND'])
        self.SCALE. append(new_dv['SCALE'])
        self.MARKER.append(new_dv['MARKER'])
        self.PARAM. append(new_dv['PARAM'])
    
    def extend(self,new_dvs):
        assert isinstance(new_dvs,DEFINITION_DV) , 'input must be of type DEFINITION_DV'
        self.KIND.  extend(new_dvs['KIND'])
        self.SCALE. extend(new_dvs['SCALE'])
        self.MARKER.extend(new_dvs['MARKER'])
        self.PARAM. extend(new_dvs['PARAM'])