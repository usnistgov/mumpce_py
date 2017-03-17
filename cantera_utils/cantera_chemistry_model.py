from state_definition import StateDefinition
from abc import ABCMeta, abstractmethod
import cantera as ct
#import mumpce

#This is added because mumpce may not be in the path and we know that mumpce exists upstairs from cantera_chemistry_model
#This line is needed for Sphinx autodoc to work. You may need to remove it yourself
import sys
sys.path.append('..')

import mumpce_py as mumpce
import numpy as np

import tqdm


class CanteraChemistryModel(mumpce.Model):
    """A class for Cantera chemistry models.
    
    :param T: The unburned gas temperature in Kelvins
    :param Patm: The pressure in atmospheres (will be converted internally to Pa)
    :param composition: The composition of the unburned gas. Can be a float array or a Cantera composition string
    :param chemistry_model: The chemistry model for the flame. Must be a chemistry model that can be used to make a Cantera phase object
    :type T: float
    :type Patm: float
    :type composition: str,ndarray(float)
    :type chemistry_model: str
    
    This is a class that implements certain common methods for the Cantera-based simulations in this package. It is a subclass of :py:class:`.Model`. It implements the following methods required by :py:class:`.Model`:
    
    * :func:`sensitivity`
    * :func:`get_parameter`
    * :func:`perturb_parameter`
    * :func:`reset_model`
    * :func:`get_model_parameter_info`
    
    You cannot instantiate a member of this class because it does not have an :func:`evaluate` method.
    
    It also includes the following methods that permit easier access to the Cantera objects that the model contains:
    
    * :func:`prepare_chemistry`
    * :func:`initialize_chemistry`
    * :func:`blank_chemistry`
    
    In addition, the method :func:`initialize_reactor` is provided as an abstract method. Subclasses of this class must define this method or else they cannot be instantiated.
    

    """
    
    __metaclass__ = ABCMeta
    
    def __init__(self,
                 T,Patm,composition,
                 chemistry_model,**kwargs):
        #Initialize the reactor's intial state
        self.initial = StateDefinition(T,Patm,composition)
        
        # Set the thermochemical and reactor models
        self.chemistry_model = chemistry_model
        
        # Create the sensitivity flag which will tell the evaluate method whether we are in a sensitivity calculation or not
        self._sens_flag = False
        
        # Initialize the Cantera mixture, thermo, and chemical model
        self.prepare_chemistry(**kwargs)
        self.tqfunc = tqdm.tqdm
        
        return
        
#    def __str__(self):
#        modelstr = str(self.initial.T) + ' K, ' + str(self.initial.P) + ' Pa ' + str(self.initial.composition)
#        return modelstr
        
    
    def prepare_chemistry(self,no_efficiencies=True,no_energy=True,no_falloff=True,**kwargs):
        """Instantiate the Cantera chemistry model and get information about the reaction model. This is called during instantiation of the model and normally would not be called at any other time.
        """
        #Call the blank_chemistry function in order to create the chemistry and simulation attributes, initialized to None
        self.blank_chemistry()
        
        #Create the Cantera gas object from the chemistry model and set its initial state
        self.initialize_chemistry()
        
        #Get the parameters that will be investigated for sensitivity analysis and find how many there are
        #self.model_parameter_info = self.get_model_parameter_info(no_efficiencies=True,no_energy=True,no_falloff=True)
        self.model_parameter_info = self.get_model_parameter_info(no_efficiencies=no_efficiencies,
                                                                  no_energy=no_energy,
                                                                  no_falloff=no_falloff)
        self.number_parameters = len(self.model_parameter_info)
        
        #Blank the chemistry so that the model can be pickled
        self.blank_chemistry()
    
    def initialize_chemistry(self):
        """Create the Cantera phase object and set its initial state
        
        This method does the following:
        
        * Checks to see if a Cantera Solution object exists that defined the thermodynamic state of the Cantera reactor, and creates that object if it does not exist
        * Sets the state of the Cantera Solution object to the state specified in self.initial
        
        """
        #If the gas object is blank, create the Cantera solution object
        if self.gas is None:
            self.gas = ct.Solution(self.chemistry_model)
        #Set the gas initial condition
        self.gas.TPX = self.initial.T, self.initial.P, self.initial.composition
        return
    
    @abstractmethod
    def initialize_reactor(self):
        """Initialize the Cantera reactor and simulation objects. This must have the following form::
        
           self.reactor = ct.SomeReactorClass(*args)
           self.simulation = ct.SomeSimulationClass(*args)
           
        This method is an abstract method and so this method must be defined by a subclass.
        
        """
        pass
    
    def blank_chemistry(self):
        """Erase the Cantera phase object, Cantera reactor object, and Cantera simulation object. 
        """
        #Set  the gas, reactor, and simulation attributes to None
        # This is used upon initializing the model object to create the attribute and before pickling or otherwise saving the object
        self.gas = None
        self.reactor = None
        self.simulation = None
    
    def load_restart(self,filename=None,solution_name=None):
        """Load a previously-saved solution from a restart file.
        
        This is a blank function intended as a placeholder for loading or saving restart files. It is defined here because :func:`sensitivity` calls these functions. The :func:`flame speed` class overrides this method with a more detailed method.
        """
        pass
    
    def save_restart(self,filename=None,solution_name=None,description=None):
        """Saves a solution to a restart file.
        
        This is a blank function intended as a placeholder for loading or saving restart files. It is defined here because :func:`sensitivity` calls these functions. The :func:`flame speed` class overrides this method with a more detailed method.
        """
        pass
    
    def ignore_restart(self):
        """Ignorse the solution found in any previously-saved restart file and solves the model from scratch
        
        This is a blank function intended as a placeholder for loading or saving restart files. It is defined here because :func:`sensitivity` calls these functions. The :func:`flame speed` class overrides this method with a more detailed method.
        """
        pass
    
    def get_parameter(self,parameter):
        """Retrieves a model parameter's value
        
        :param parameter_id: The parameter identifier. 
        :type parameter_id: int
        :returns: parameter_value
        :rtype: float
        """        
        param_info = self.model_parameter_info[parameter]
        reaction_number = param_info['reaction_number']
        parameter_value = self.gas.multiplier(reaction_number)
        return parameter_value
    
    def perturb_parameter(self,parameter,factor):
        """Perturbs a model parameter's value by a specified amount.
        
        :param parameter_id: The parameter identifier. 
        :type parameter_id: int
        :param new_value: The amount to change the parameters value.
        :type new_value: float
        """
        param_info = self.model_parameter_info[parameter]
        reaction_number = param_info['reaction_number']
        self.gas.set_multiplier(factor,reaction_number)
        #print param_info['parameter_name']
        return 
    
    def reset_model(self):
        """Reset all model parameters to their original values"""
        if self.gas is None:
            self.initialize_chemistry()
        self.gas.set_multiplier(1.0)
        return
    
    def get_reaction_info(self,reaction_number,reaction):
        """Gets information about a particular reaction within a Cantera model and find which parameters it has that might be active. This is called during instantiation of the model and normally would not be called at any other time.
        
        :param reaction_number: The number of the reaction within the Cantera model
        :param reaction: The reaction object corresponding to that number
        :type reaction_number: int
        :type reaction: Cantera reaction object
        :returns: reaction_info, a list of dictionaries that describe the parameters available in the model
        :rtype: list of dicts
        
        """
        
        reaction_name = self.gas.reaction_equations([reaction_number])[0]
        
        #All reactions have an A that could be active
        reaction_info  = [{'reaction_number':reaction_number,'parameter_type':'A_factor','parameter_name':reaction_name}]
        
        #Default assumption is that a reaction has no third body efficiencies or falloff behavior
        HasThirdBody = False #Default, no efficiencies
        HasFalloff = False #Default, no falloff behavior
        
        if reaction.reaction_type == 4:
            HasThirdBody = True
            HasFalloff = True
        if reaction.reaction_type == 8:
            HasThirdbody = True
            HasFalloff = True
        if reaction.reaction_type ==2:
            HasThirdBody = True
        
        if HasFalloff:
            rate = reaction.high_rate
            reaction_info  = [{'reaction_number':reaction_number,'parameter_type':'High_pressure_A','parameter_name':reaction_name}]
            #Do not consider activation energies very close to zero
            if abs(rate.activation_energy) > 0.1:
                reaction_info += [{'reaction_number':reaction_number,'parameter_type':'High_pressure_E','parameter_name':reaction_name}]
            rate = reaction.low_rate
            reaction_info += [{'reaction_number':reaction_number,'parameter_type':'Low_pressure_A','parameter_name':reaction_name}]
            if abs(rate.activation_energy) > 0.1:
                reaction_info += [{'reaction_number':reaction_number,'parameter_type':'Low_pressure_E','parameter_name':reaction_name}]
        else:
            rate = reaction.rate
            reaction_info  = [{'reaction_number':reaction_number,'parameter_type':'A_factor','parameter_name':reaction_name}]
            if abs(rate.activation_energy) > 0.1:
                reaction_info += [{'reaction_number':reaction_number,'parameter_type':'Energy','parameter_name':reaction_name}]

        if HasThirdBody:
            num_efficiencies = len(reaction.efficiencies)
            for species_name in reaction.efficiencies:
                if reaction.efficiency(species_name) > 0:
                    reaction_info += [{'reaction_number':reaction_number,'parameter_type':'Efficiency','species':species_name,
                                       'parameter_name':reaction_name}]
        return reaction_info
    
    def get_model_parameter_info(self,no_efficiencies=False,no_energy=False,no_falloff=False):
        """Gets the list of available parameters for this model
        
        :param no_efficiencies: If True, then do not consider the third-body efficiencies as active parameters
        :param no_energies: If True, then do not consider activation energies as active parameters
        :param no_falloff: If True, then do not consider high- and low-pressure limits as active parameters
        :returns: model_parameter_info
        :rtype: list of dicts
        """
        #Initialize the Cantera model
        model = ct.Solution(self.chemistry_model)
        #Initialize the model parameter info lists
        model_parameter_info_full = []
        model_parameter_info = []
        
        #Get the list of possibly-active model parameters
        for reaction_num in range(model.n_reactions):
            reaction = model.reaction(reaction_num)
            reac_info = self.get_reaction_info(reaction_num,reaction)
            model_parameter_info_full += reac_info
        for param_info in model_parameter_info_full:
            include_this_parameter = True
            if no_efficiencies:
                if param_info['parameter_type'] == 'Efficiency':
                    include_this_parameter = False
            #else:
            #    if param_info[1] == 'Efficiency':
            #        if not(param_info[2] > 0):
            #            include_this_parameter = False
            if no_energy:
                if param_info['parameter_type'] == 'Energy':
                    include_this_parameter = False
                if param_info['parameter_type'] == 'Low_pressure_E':
                    include_this_parameter = False
                if param_info['parameter_type'] == 'High_pressure_E':
                    include_this_parameter = False
            if no_falloff:
                if param_info['parameter_type'] == 'Low_pressure_A':
                    include_this_parameter = False
                if param_info['parameter_type'] == 'Low_pressure_E':
                    include_this_parameter = False
                #if param_info['parameter_type'] == 'High_pressure_A':
                #    include_this_parameter = False
                #if param_info['parameter_type'] == 'High_pressure_E':
                #    include_this_parameter = False
            if include_this_parameter:
                model_parameter_info += [param_info]
            #print param_info
        return model_parameter_info
    
    def sensitivity(self,perturbation,parameter_list,logfile):
        """Evaluates the sensitivity of the model value with respect to the model parameters
        
        :param perturbation: The amount to perturb each parameter during the sensitivity analysis
        :param parameter_list: The list of parameters to perturb. This will be a list of parameter identifiers, which are usually ints or strs.
        :param logfile: The logging file that will contain the sensitivity calculation output.
        :type perturbation: float
        :type parameter_list: array_like
        :type logfile: str
        :returns: model_value,sensitivity_vector
        :rtype: float,ndarray
        """
        #Intialize the sensitivity vector
        sensitivity_vector = np.zeros(len(parameter_list))
        
        #Evaluate the model once and save the result in a restart file
        value = self.evaluate()
        self.save_restart()
        #print("Value = {: 10.5e}".format(value))
        logfile.write("Value = {: 10.5e}".format(value))
        
        pos_mult = 1 + perturbation
        neg_mult = 1/pos_mult
        
        #print pos_mult
        #print neg_mult
        
        logfile.write('Rxn  Value+       Value-           Sensitivi   Reaction Name')
        
        #pbar = tqdm.tqdm(total=len(parameter_list))
        
        self._sens_flag = True
        for (param_number,param_id) in enumerate(self.tqfunc(parameter_list,desc=logfile.name)):
        #for (param_number,param_id) in enumerate(parameter_list):
            #pbar.update(1)
            

            mult_base = self.get_parameter(param_id)
            #print mult_base
            pos_pert = pos_mult*mult_base
            #print pos_pert
            self.perturb_parameter(param_id,pos_pert)
            #print "going into ignition delay problem"
            valuep = self.evaluate()
            #for parmid in parameter_list:
            #    print parmid,self.get_parameter(parmid)
            neg_pert = neg_mult*mult_base
            #print neg_pert
            self.perturb_parameter(param_id,neg_pert)
            self.load_restart()
            valuem = self.evaluate()
            #for parmid in parameter_list:
            #    print parmid,self.get_parameter(parmid)            
            self.perturb_parameter(param_id,mult_base)
            
            #sensitivity = (delayp - delaym) / (2.0 * perturbation * delay)
            
            sensitivity_vector[param_number] = (valuep - valuem) / (2.0 * perturbation * value)
            
            #print (delayp - delaym)
            #print delay
            #print (delayp - delaym) / delay
            #print (delayp - delaym) / (2.0 * perturbation * delay)
            #print sensitivity
            #print perturbation
            #print sensitivity_vector[param_number]
            #print('{: 4d} {: 10.5e}  {: 10.5e}  {: 10.4e}  {}'.format(param_id,
            #                                                      valuep,valuem,sensitivity_vector[param_number],
            #                                                      self.gas.reaction_equations([param_id])[0])
            #       )
            logfile.write('{: 4d} {: 10.5e}  {: 10.5e}  {: 10.4e}  {}\n'.format(param_id,
                                                                  valuep,valuem,sensitivity_vector[param_number],
                                                                  self.gas.reaction_equations([param_id])[0])
                   )
        #value = math.log(value/1.0e-6)
        self._sens_flag = False
        return value, sensitivity_vector
    
    def print_sens(self,sensitivity_vector,print_params=None):
        """Sorts the parameters by sensitivity coefficient and prints them. Specify max_number to print only max_number parameters
        """
        
        
        if print_params is None:
            print_params = self.number_parameters
        
        for print_param in print_params:
            print('{: 4d} {: 10.4e}  {}'.format(print_param,
                                                sensitivity_vector[print_param],
                                                self.gas.reaction_equations([print_param])[0]
                                               )
                 )
        
        return