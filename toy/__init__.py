import sys
sys.path.append('..')

import mumpce_py as mumpce
#import mumpce
import numpy as np
import pandas as pd

zeros         = np.array([[0.2 ,-0.6,0.5 ,0.8 ,0.5 ]]).transpose()
zeros_app     = np.array([[0.0, 0.0, 0.0]]).transpose()

first_order = np.array([[1.00, 0.02, 1.00, 0.40, 0.03 ],
                        [0.60, 0.60, 0.40, 1.00, 0.10 ],
                        [0.40, 1.00, 0.60, 0.01, 1.00 ],
                        [0.10, 0.03, 0.10, 0.02, 0.02 ],
                        [0.03, 0.10, 0.02, 0.60, 0.60 ],
                        [0.01, 0.01, 0.03, 0.03, 0.01 ],
                        [0.02, 0.40, 0.01, 0.10, 0.40 ]]
                      )

first_order_app = np.array([[1.00, 0.01, 0.01 ],
                            [0.20, 0.20, 0.01 ],
                            [0.01, 1.00, 0.01 ],
                            [0.01, 0.20, 0.01 ],
                            [0.01, 0.01, 0.20 ],
                            [0.01, 0.01, 1.00 ],
                            [0.01, 0.01, 0.20 ]]
                          )

values        = np.array([[0.4 ,0.25,0.2 ,0.7 ,0.4 ]]).transpose()
uncertainties = np.array([[0.05,0.05,0.05,0.08,0.8 ]]).transpose()
parameter_uncertainties = np.exp(np.ones(7))

class toy_model(mumpce.Model):
    """An example of a generic model for use with the MUMPCE program.
    
    This is a toy model where the response is given by :math:`\ln y = ax + z`, where :math:`x` is the vector of parameters that can be perturbed, :math:`a` is some sensitivity vector and :math:`z` is some zero-order term. All of the parameters are nominally set to zero, with +/- 1 representing the upper and lower bounds of the prior uncertainty for these parameters.
    
    The model has a matrix of five possible hard-wired :math:`z` and :math:`a` values. The experiment number corresponds to which pair is being used when the model is being instantiated.
    
    :param experiment_number: Which experiment (of five possible) this model represents
    :type experiment_number: int
    :returns: A MUM-PCE model object representing this measurement
    """
    def __init__(self,experiment_number,loglevel=True):
        self.experiment_number = experiment_number
        
        self.number_parameters = 7
        
        self.parameter_vector = np.zeros((self.number_parameters,1))
        
        self.model_parameter_info = self.get_model_parameter_info()
        
        self.loglevel = loglevel
        
    def __str__(self):
        """Returns the experiment number, which is the only thing interesting about this model
        """
        return str(self.experiment_number)
    
    def evaluate(self):
        """Run the model once and return a single value
        
        :returns: model_value
        """
        #This is the set of Z terms. There are seven possible terms 
        zero_term = zeros
        a_term = first_order
        #first_order = np.array([[1.00, 0.02, 1.00, 0.40, 0.03 ],
        #                        [0.60, 0.60, 0.40, 1.00, 0.10 ],
        #                        [0.40, 1.00, 0.60, 0.01, 1.00 ],
        #                        [0.10, 0.03, 0.10, 0.02, 0.02 ],
        #                        [0.03, 0.10, 0.02, 0.60, 0.60 ],
        #                        [0.01, 0.01, 0.03, 0.03, 0.01 ],
        #                        [0.02, 0.40, 0.01, 0.10, 0.40 ]]
        #                      )
        
        #print zero_term
        #print np.dot(first_order.T,self.parameter_vector)
        
        value_vector = zero_term + np.dot(a_term.T,self.parameter_vector)
        #print value_vector
        value = value_vector[self.experiment_number]
        #print value
        return np.exp(value[0])
    
    def sensitivity(self,perturbation=1.0e-3,parameter_list=None,logfile='file'):
        """Evaluate the sensitivity of the model value with respect to the model parameters
        
        :param perturbation: The amount to perturb each parameter during the sensitivity analysis
        :param parameter_list: The list of parameters to perturb. Usually this will be a list of parameter identifiers, which are usually ints or strs.
        :param logfile: The logging file that will contain the sensitivity calculation output.
        :type perturbation: float
        :type parameter_list: iterable
        :type logfile: str
        :returns: model_value,sensitivity_vector
        """
        
        if parameter_list is None:
            parameter_list = range(self.number_parameters)
        
        sensitivity_vector = np.zeros(len(parameter_list))
        
        value = self.evaluate()
        if self.loglevel: print ("Value = ", value)
        
        pos_mult = 1 + perturbation
        neg_mult = 1/pos_mult
        
        if self.loglevel: print ('Param  Value+       Value-           Sensitivity')
        
        for (param_number,param_id) in enumerate(parameter_list):
            #Get the base parameter value
            mult_base = self.get_parameter(param_id)
            
            #Calculate the positive perturbation, perturb the model, and evaluate the model
            pos_pert = pos_mult*mult_base
            self.perturb_parameter(param_id,pos_pert)
            valuep = self.evaluate()
            
            #Calculate the negative perturbation, perturb the model, and evaluate the model
            neg_pert = neg_mult*mult_base
            self.perturb_parameter(param_id,neg_pert)
            valuem = self.evaluate()            
            self.perturb_parameter(param_id,mult_base)

            
            sensitivity_vector[param_number] = (valuep - valuem) / (2.0 * perturbation * value)
            
            if self.loglevel: 
                print(
                    '{: 6d} {: 10.5e}  {: 10.5e}  {: 10.4e}'.format(param_id,valuep,valuem,sensitivity_vector[param_number])
                     )
        
        return value,sensitivity_vector
    
    def get_parameter(self,parameter_id):
        """Retrieve a model parameter's value
        
        :param parameter_id: The parameter number whose value to retrieve
        :type parameter_id: int
        :returns: parameter_value
        :rtype: float
        """
        parameter_value = np.exp(self.parameter_vector[parameter_id])
        return parameter_value
    
    def perturb_parameter(self,parameter_id,factor):
        """Perturb a model parameter's value by a specified amount. 
        
        :param parameter_id: The parameter identifier whose value to modify
        :param factor: The amount to change the parameter's value.
        :type factor: float
        """
        self.parameter_vector[parameter_id] = np.log(factor)
        return
    
    def reset_model(self):
        """Reset all model parameters to their original values"""
        self.parameter_vector = np.zeros((7,1))
        return
    
    def get_model_parameter_info(self):
        """Get information about the parameters, which will go up to the hosting measurement. This is called during instantiation of the model and normally would not be called at any other time. 

        :returns: model_parameter_info
        :rtype: list of dicts
        """
        model_parameter_info = []
        for parameter_number in range(self.number_parameters):
            param_info = [{'parameter_number':parameter_number,'parameter_name':'Parameter ' + str(parameter_number+1)}]
            model_parameter_info += param_info
        model_parameter_info = np.array(model_parameter_info)
        return model_parameter_info

class toy_app(toy_model):
    """An example of a generic model for use with the MUMPCE program.
    
    This is identical to :func:`toy.toy_model`, except that it is intended as an targeted application to highlight the experimental design functions. There are three possible hard-wired :math:`z` and :math:`a` values, which are different from those available in :func:`toy.toy_model`
    
    :param experiment_number: Which targeted application (of three possible) this model represents
    :type experiment_number: int
    :returns: A MUM-PCE model object representing this experiment
    """
    def evaluate(self):
        
        zero_term = zeros_app
        a_term = first_order_app
        #first_order = np.array([[1.00, 0.01, 0.01 ],
        #                        [0.20, 0.20, 0.01 ],
        #                        [0.01, 1.00, 0.01 ],
        #                        [0.01, 0.20, 0.01 ],
        #                        [0.01, 0.01, 0.20 ],
        #                        [0.01, 0.01, 1.00 ],
        #                        [0.01, 0.01, 0.20 ]]
        #                      )
        
        #print zero_term
        #print np.dot(first_order.T,self.parameter_vector)
        
        value_vector = zero_term + np.dot(a_term.T,self.parameter_vector)
        #print value_vector
        value = value_vector[self.experiment_number]
        #print value
        return np.exp(value[0])

def toy_initialize(filename,model):
    """The initialization function for the toy model. This function shows how to instantiate a MUMPCE model object and then associate it with a MUMPCE measurement object. It will read an experimental database from an Excel file and use this information to build the measurement lists.
    
    :param filename: The file that contains the experimental database.
    :type filename: str
    :returns: measurement_list, a list of MUMPCE measurement objects
    :rtype: list
    """
    
    #Create the blank measurement list
    measurement_list = []
    
    #Read the Excel data for this project
    df = pd.read_excel(filename)
    
    #Get the list of column names
    df_columns = df.columns.values
    
    #Group the experiments by unique ID so that we can iterate over them
    gb = df.groupby('Name')
    
    for name,this_experiment in gb:
        exp_name = this_experiment.Name.values[0]
        exp_num = this_experiment.Number.values[0]
        
        #Use a different model depending on whether this is an Experiment or an Application
        if exp_name.startswith('Exp'):
            this_model = toy_model(exp_num)
        else:
            this_model = toy_app(exp_num)
        
        this_value = None
        this_unc = None
        
        if 'Value' in df_columns:
            this_value = this_experiment.Value.values[0]
            this_unc = this_experiment.Uncertainty.values[0]
        
        #Call mumpce.measurement to create the measurement object
        meas = mumpce.Measurement(name=exp_name,
                                  model=this_model,
                                  value=this_value,
                                  uncertainty=this_unc,
                                  response_perturbation=0.0000001,
                                  response_type='log'
                                 )
        #Add the measurement object to the list
        measurement_list += [meas]
    return measurement_list