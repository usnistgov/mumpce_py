import cantera as ct
#from cantera import one_atm

class StateDefinition(object):
    """A top level class for a thermodynamic state definition of an ideal gas
    
    :param T: The unburned gas temperature in Kelvins
    :param Patm: The pressure in atmospheres (will be converted internally to Pa)
    :param composition: The composition of the unburned gas. Can be a float array or a Cantera composition string
    :type T: float
    :type Patm: float
    :type composition: str,ndarray(float)
    
    """
    def __init__(self,
                T,Patm,composition):
        self.T = T
        self.P = Patm * ct.one_atm#ct.one_atm
        self.composition = composition
        return
    
#    def __str__(self):
#        
