from cantera_chemistry_model import CanteraChemistryModel
import numpy as np
import cantera as ct
import mumpce_py as mumpce
from mumpce_py.response_surface import ResponseSurface

class RxnMeasurement(mumpce.Measurement):
    """A special class for optimizing reaction rate measurements.
    
    Since there is an analytic expression for the sensitivities, creating the response surfaces is more straigtforward than it would be for a normal experiment. The call signature for this measurement is the same as for :py:class:`mumpce.Measurement`.
    """
    
    def make_response(self):
        """Generates a sensitivity_analysis_based response surface for this measurement
        """
        #zero_term = self.evaluate
        self.model.reset_model()
        logfile_name = self.name + '_resp_log.out'
        response_logfile = open(logfile_name,'w')
        
        sensitivity_args = (self.response_sensitivity,
                           self.active_parameters,
                           response_logfile)
        
        zero_term, sens_zero = self.model.sensitivity(*sensitivity_args)
        
        if self.response_type == 'log':
            zero_term = np.log(zero_term)

        
        number_params = len(self.active_parameters)
        a_terms = np.zeros((number_params,1))
        b_terms = np.zeros((number_params,number_params))
        d_terms = np.zeros_like(b_terms)
        
        for parameter_number,(parameter,sens_val,uncert) in enumerate(zip(self.active_parameters,
                                                                              sens_zero,
                                                                              self.parameter_uncertainties)):
            parameter_type = self.model.model_parameter_info[parameter]['parameter_type']
            if 'A' in parameter_type:
                a_terms[parameter_number] = sens_val * np.log(uncert)
            if 'E' in parameter_type:
                a_terms[parameter_number] = sens_val * (uncert - 1)
        
        self.response = ResponseSurface(zero_term=zero_term,
                                         a_terms=a_terms,
                                         b_terms=b_terms,
                                         d_terms=d_terms,
                                         active_parameters=self.active_parameters)
        return

class ReactionRateBase(CanteraChemistryModel):
    """The base class for optimizing reaction rates.
    
    This class defines some shared attributes across all ReactionRate type classes.
    
    :param T: The temperature at which to evaluate the rate constant in Kelvins
    :param Patm: The pressure in atmospheres (will be converted internally to Pa)
    :param composition: The composition for which the rate constant will be evaluated. This is important for reactions with third-body collision efficiencies. Can be a float array or a Cantera composition string
    :type T: float
    :type Patm: float
    :type composition: str,ndarray(float)
    
    """
    def __init__(self,
                 T,Patm,composition,
                 chemistry_model,
                 **kwargs):
        
        #Generic interface to CanteraChemistryModel
        super(ReactionRateBase,self).__init__(T,Patm,composition,chemistry_model,**kwargs)
        
        self.parameter_list = []
        self.numerator_list = []
        self.denominator_list = []
        
        return
    
    def sensitivity(self,perturbation,parameter_list,logfile):
        sensitivity_list = []
        
        value = self.evaluate()
        logfile.write("Value = {: 10.5e}\n".format(value))
        
        pos_mult = 1 + perturbation
        neg_mult = 1/pos_mult
        
        for (param_number,param_id) in enumerate(self.tqfunc(parameter_list,desc=logfile.name)):
            
            #Default value for sensitivity is 0, because most parameters are rate parameters for other reactions
            sensitivity = 0.0
            
            #Check to see if this parameter is part of this reaction's rate expression
            if param_id in self.parameter_list:
                
                #Get what kind of parameter this is
                param_info = self.model_parameter_info[param_id]
                parameter_type = param_info['parameter_type']
                
                #Simple Arrhenius expression has a very simple form for the sensitivity
                if 'pressure' not in parameter_type:
                #if False:
                    if 'A' in parameter_type:
                        sensitivity = 1.0
                    if 'E' in parameter_type:
                        E = self.get_parameter(param_id)
                        sensitivity = -E/(ct.gas_constant*self.initial.T)
                else: #More complex reaction -> do brute sensitivity
                    mult_base = self.get_parameter(param_id)
                    
                    #positive perturbation
                    pos_pert = pos_mult*mult_base
                    self.perturb_parameter(param_id,pos_pert)
                    valuep = self.evaluate()
                    
                    #negative perturbation
                    neg_pert = neg_mult*mult_base
                    self.perturb_parameter(param_id,neg_pert)
                    valuem = self.evaluate()
                    
                    #restore original value
                    self.perturb_parameter(param_id,mult_base)
                    
                    sensitivity = (valuep - valuem) / (2.0 * perturbation * value)
                if param_id in self.denominator_list:
                    sensitivity = sensitivity * -1.0
            
            sensitivity_list += [sensitivity]
        
        sensitivity_vector = np.array(sensitivity_list)
        
        return value,sensitivity_vector
    
    def initialize_reactor():
        pass
        
        
class ReactionRateAtCondition(ReactionRateBase):
    """A class that will specify a rate constant at a particular temperature and pressure.
    
    :param T: The temperature at which to evaluate the rate constant in Kelvins
    :param Patm: The pressure in atmospheres (will be converted internally to Pa)
    :param composition: The composition for which the rate constant will be evaluated. This is important for reactions with third-body collision efficiencies. Can be a float array or a Cantera composition string
    :key reaction_number: The reaction number from the Cantera model which will be used
    :type T: float
    :type Patm: float
    :type composition: str,ndarray(float)
    :type reaction_number: int
    
    """
    
    def __init__(self,
                 T,Patm,composition,
                 chemistry_model,
                 reaction_number=None,**kwargs):
        
        super(ReactionRateAtCondition,self).__init__(T,Patm,composition,chemistry_model,**kwargs)
        
        self.rxn_num = reaction_number
        self.rxn_name = ''
        
        for (param_number,param_info) in enumerate(self.tqfunc(self.model_parameter_info)):
            if param_info['reaction_number'] == reaction_number:
                self.parameter_list += [param_number]
                self.rxn_name = param_info['parameter_name']
                
        return
    
    def __str__(self):
        
        str_args = (self.rxn_name,
                    self.initial.T,                    
                    self.initial.P/1.0e3,
                   )
        
        modelstr = 'Reaction {}: {:8.0f} K, {:5.2f} kPa'.format(*str_args)
        return modelstr
    
    def evaluate(self):
        self.initialize_chemistry()
        
        return self.gas.forward_rate_constants[self.rxn_num]
        
        #return
    

class ReactionRateRatioAtCondition(ReactionRateBase):
    """A class that will specify a ratio of two rate constants at a particular temperature and pressure.
    
    :param T: The temperature at which to evaluate the rate constant in Kelvins
    :param Patm: The pressure in atmospheres (will be converted internally to Pa)
    :param composition: The composition for which the rate constant will be evaluated. This is important for reactions with third-body collision efficiencies. Can be a float array or a Cantera composition string
    :key reaction_numerator: The reaction number from the Cantera model that specifies the numerator reaction
    :key reaction_denominator: The reaction number from the Cantera model that specifies the denominator reaction
    :type T: float
    :type Patm: float
    :type composition: str,ndarray(float)
    :type reaction_numerator: int
    :type reaction_denominator: int
    """
    
    def __init__(self,
                 T,Patm,composition,
                 chemistry_model,
                 reaction_numerator=None,reaction_denominator=None,**kwargs):
        
        super(ReactionRateRatioAtCondition,self).__init__(T,Patm,composition,chemistry_model,**kwargs)
        
        self.rxn_num = reaction_numerator
        self.rxn_den = reaction_denominator
        
        for (param_number,param_info) in enumerate(self.tqfunc(self.model_parameter_info)):
            if param_info['reaction_number'] == reaction_numerator:
                self.numerator_list += [param_number]
                self.numer_name = param_info['parameter_name']
        
        for (param_number,param_info) in enumerate(self.tqfunc(self.model_parameter_info)):
            if param_info['reaction_number'] == reaction_denominator:
                self.denominator_list += [param_number]
                self.denom_name = param_info['parameter_name']
        self.parameter_list = self.numerator_list + self.denominator_list
        
        return
    
    def __str__(self):
        
        str_args = (self.numer_name,self.denom_name,
                    self.initial.T,                    
                    self.initial.P/1.0e3,
                   )
        
        modelstr = 'Rate ratio {}/{}: {:8.0f} K, {:5.2f} kPa'.format(*str_args)
        return modelstr
    
    def evaluate(self):
        self.initialize_chemistry()
        
        numer = self.gas.forward_rate_constants[self.rxn_num]
        denom = self.gas.forward_rate_constants[self.rxn_den]
        
        return numer/denom   
    
class ReactionARatio(ReactionRateBase):
    """A class that will fix the ratio of A-factors for two reactions
    
    :param T: The temperature at which to evaluate the rate constant in Kelvins
    :param Patm: The pressure in atmospheres (will be converted internally to Pa)
    :param composition: The composition for which the rate constant will be evaluated. This is important for reactions with third-body collision efficiencies. Can be a float array or a Cantera composition string
    :key reaction_numerator: The reaction number from the Cantera model that specifies the numerator reaction
    :key reaction_denominator: The reaction number from the Cantera model that specifies the denominator reaction
    :type T: float
    :type Patm: float
    :type composition: str,ndarray(float)
    :type reaction_numerator: int
    :type reaction_denominator: int
    
    """
    def __init__(self,
                 T,Patm,composition,
                 chemistry_model,
                 reaction_numerator=None,reaction_denominator=None,**kwargs):
        
        super(ReactionARatio,self).__init__(T,Patm,composition,chemistry_model,**kwargs)
        
        self.rxn_num = reaction_numerator
        self.rxn_den = reaction_denominator
        
        print(reaction_numerator,reaction_denominator)
        
        #self.rxn = 
        
        self.numerator_num = []
        self.denominator_list = []
        
        for (param_number,param_info) in enumerate(self.tqfunc(self.model_parameter_info)):
            if 'A' in param_info['parameter_type']:
                if param_info['reaction_number'] == reaction_numerator:
                    self.numer_num = param_number
                    self.numer_name = param_info['parameter_name']
                if param_info['reaction_number'] == reaction_denominator:
                    self.denom_num = param_number
                    self.denom_name = param_info['parameter_name']
        
        self.parameter_list = [self.numer_num] + [self.denom_num]
        self.denominator_list = [self.denom_num]
        
        return
    
    def __str__(self):
        
        str_args = (self.numer_name,self.denom_name,
                   )
        
        modelstr = 'A factor ratio ({})/({})'.format(*str_args)
        return modelstr
    
    def evaluate(self):
        self.initialize_chemistry()
        
        numer = self.get_parameter(self.numer_num)
        denom = self.get_parameter(self.denom_num)
        
        return numer/denom

class ReactionEDiff(ReactionRateBase):
    """A class that will fix the difference between activation energies for two reactions
    
    :param T: The temperature at which to evaluate the rate constant in Kelvins
    :param Patm: The pressure in atmospheres (will be converted internally to Pa)
    :param composition: The composition for which the rate constant will be evaluated. This is important for reactions with third-body collision efficiencies. Can be a float array or a Cantera composition string
    :key reaction_numerator: The reaction number from the Cantera model that specifies the numerator reaction
    :key reaction_denominator: The reaction number from the Cantera model that specifies the denominator reaction
    :type T: float
    :type Patm: float
    :type composition: str,ndarray(float)
    :type reaction_numerator: int
    :type reaction_denominator: int
    
    """
    def __init__(self,
                 T,Patm,composition,
                 chemistry_model,
                 reaction_numerator=None,reaction_denominator=None,**kwargs):
        
        T = 1#/ct.gas_constant
        
        super(ReactionEDiff,self).__init__(T,Patm,composition,chemistry_model,**kwargs)
        
        self.rxn_num = reaction_numerator
        self.rxn_den = reaction_denominator
        
        print(reaction_numerator,reaction_denominator)
        
        #self.rxn = 
        
        self.numerator_num = []
        self.denominator_list = []
        
        for (param_number,param_info) in enumerate(self.tqfunc(self.model_parameter_info)):
            if 'E' in param_info['parameter_type']:
                if param_info['reaction_number'] == reaction_numerator:
                    self.numer_num = param_number
                    self.numer_name = param_info['parameter_name']
                if param_info['reaction_number'] == reaction_denominator:
                    self.denom_num = param_number
                    self.denom_name = param_info['parameter_name']

        
        self.parameter_list = [self.numer_num] + [self.denom_num]
        self.denominator_list = [self.denom_num]
        
        return
    
    def __str__(self):
        
        str_args = (self.numer_name,self.denom_name,
                   )
        
        modelstr = 'Activation energy difference ({})/({})'.format(*str_args)
        return modelstr
    
    def evaluate(self):
        self.initialize_chemistry()
        
        numer = self.get_parameter(self.numer_num)
        denom = self.get_parameter(self.denom_num)
        
        return np.exp(numer - denom)
