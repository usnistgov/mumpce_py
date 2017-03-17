from cantera_chemistry_model import CanteraChemistryModel
import numpy as np
import cantera as ct

class ShockTube(CanteraChemistryModel):
    """A class for defining shock tube simulations.
    
    This is a class that will create a simulation of a shock tube. The shock tube simulations are subclasses of this class. It is a subclass of :func:`cantera_chemistry_model`, which is in turn a subclass of :func:`model`. You cannot instantiate a member of this class because it does not have an :func:`evaluate` method.
    
    This class implements the :func:`initialize_reactor` method required by :func:`cantera_chemistry_model`. It also redefines :func:`reset_model`.
    
    :param T: The unburned gas temperature in Kelvins
    :param Patm: The pressure in atmospheres (will be converted internally to Pa)
    :param composition: The composition of the unburned gas. Can be a float array or a Cantera composition string
    :param reactor_model: The Cantera reactor model used to simulate the shock tube. Typically IdealGasReactor
    :param chemistry_model: The chemistry model for the shock tube. Must be a filename that contains a Cantera chemistry model that can be used to make a Cantera phase object
    :param loglevel: The level of external logging that the model should do. If None, does not produce any logging in the ignition delay solver. Otherwise, produces detailed output
    :type T: float
    :type Patm: float
    :type composition: array-like or str
    :type reactor_model: Cantera reactor model
    :type chemistry_model: str
    """    
    def __init__(self,
                 T,Patm,composition,
                 reactor_model,chemistry_model,loglevel=None,**kwargs):
        
        super(ShockTube,self).__init__(T,Patm,composition,chemistry_model,**kwargs)
        self.reactor_model = reactor_model #: The Cantera :py:class:`cantera.IdealGasReactor` class that will be instantiated by :py:func:`initalize_reactor`
        self.loglevel = loglevel
        
        return
    
    def initialize_reactor(self):
        """Initialize the shock tube. It consists of a Cantera :py:class:`cantera.IdealGasReactor` object defined in self.reactor_model and a Cantera reactor network containing only that reactor.
        """        
        self.reactor = self.reactor_model(self.gas)
        self.simulation = ct.ReactorNet([self.reactor])
        
        self.simulation.atol = 1.0e-13
        self.simulation.rtol = 1.0e-4
        
        
        return
    
    def get_parameter_thisisthecomplicatedonethatdoesntwork(self,parameter):
        param_info = self.model_parameter_info[parameter]
        reaction_number = param_info['reaction_number']
        parameter_type = param_info['parameter_type']
        
        reaction = self.gas.reaction(reaction_number)
        rate = None
        efficiencies = None
        
        if parameter_type == 'A_factor':
            rate = reaction.rate
            parameter_value = rate.pre_exponential_factor
        if parameter_type == 'Energy':
            rate = reaction.rate
            parameter_value = rate.activation_energy
        if parameter_type == 'High_pressure_A':
            rate = reaction.high_rate
            parameter_value = rate.pre_exponential_factor
        if parameter_type == 'High_pressure_E':
            rate = reaction.high_rate
            parameter_value = rate.activation_energy
        if parameter_type == 'Low_pressure_A':
            rate = reaction.low_rate
            parameter_value = rate.pre_exponential_factor
        if parameter_type == 'Low_pressure_E':
            rate = reaction.low_rate
            parameter_value = rate.activation_energy
        if parameter_type == 'Efficiency':
            efficiencies = reaction.efficiencies
            species = param_info['species']
            parameter_value = efficiencies[species]
        return parameter_value
    
    def perturb_parameter_thisisthecomplicatedonethatdoesntwork(self,parameter,new_value):
        param_info = self.model_parameter_info[parameter]
        reaction_number = param_info['reaction_number']
        parameter_type = param_info['parameter_type']
        
        print parameter
        print param_info
        
        reaction = self.gas.reaction(reaction_number)
        print reaction.rate
        rate = None
        efficiencies = None
        reactants = reaction.reactant_string
        products = reaction.product_string
        
        rtype = reaction.reaction_type
        
        pressurestring = 'pressure'
        HasFalloff = False
        
        if pressurestring in parameter_type:
            HasFalloff = True
        if HasFalloff:
            highrate = reaction.high_rate
            lowrate = reaction.low_rate
            if rtype == 4:
                newreaction = ct.FalloffReaction(reactants=reactants,products=products)
            if rtype == 8:
                newreaction = ct.ChemicallyActivatedReaction(reactants=reactants,products=products)
        else:
            rate = reaction.rate
            if rtype == 2:
                newreaction = ct.ThreeBodyReaction(reactants=reactants,products=products)
            else:
                newreaction = ct.ElementaryReaction(reactants=reactants,products=products)
        
            
        
        if parameter_type == 'A_factor':
            old_A = rate.pre_exponential_factor
            old_b = rate.temperature_exponent
            old_E = rate.activation_energy
            new_A = new_value
            newreaction.rate = ct.Arrhenius(A=new_A,b=old_b,E=old_E)
        if parameter_type == 'Energy':
            old_A = rate.pre_exponential_factor
            old_b = rate.temperature_exponent
            old_E = rate.activation_energy
            new_E = new_value
            newreaction.rate = ct.Arrhenius(A=old_A,b=old_b,E=new_E)
        if parameter_type == 'High_pressure_A':
            rate = highrate
            old_A = rate.pre_exponential_factor
            old_b = rate.temperature_exponent
            old_E = rate.activation_energy
            new_A = new_value
            newreaction.high_rate = ct.Arrhenius(A=new_A,b=old_b,E=old_E)
        if parameter_type == 'High_pressure_E':
            rate = highrate
            old_A = rate.pre_exponential_factor
            old_b = rate.temperature_exponent
            old_E = rate.activation_energy
            new_E = new_value
            newreaction.high_rate = ct.Arrhenius(A=old_A,b=old_b,E=new_E)
        if parameter_type == 'Low_pressure_A':
            rate = lowrate
            old_A = rate.pre_exponential_factor
            old_b = rate.temperature_exponent
            old_E = rate.activation_energy
            new_A = new_value
            newreaction.low_rate = ct.Arrhenius(A=new_A,b=old_b,E=old_E)
        if parameter_type == 'Low_pressure_E':
            rate = lowrate
            old_A = rate.pre_exponential_factor
            old_b = rate.temperature_exponent
            old_E = rate.activation_energy
            new_E = new_value
            newreaction.low_rate = ct.Arrhenius(A=old_A,b=old_b,E=new_E)
            
        self.gas.modify_reaction(reaction_number,newreaction)
            
        if parameter_type == 'Efficiency':
            efficiencies = copy.deepcopy(reaction.efficiencies)
            species = param_info['species']
            efficiencies[species] = new_value
            reaction.efficiencies = efficiencies
        return 
    
    def reset_model(self):
        """Reset all model parameters to their original values
        
        This version works by erasing the chemistry and re-initializing it from the CTI file.
        """
        self.blank_chemistry()
        self.initialize_chemistry()
        return
