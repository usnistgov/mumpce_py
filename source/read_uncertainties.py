import cantera as ct
import numpy as np

def read_uncertainties(uncertainty_file=None,mumpce_cantera_model=None):
    #Read the uncertainty file
    with open(uncertainty_file,'r') as uncert_file:
        uncerts = np.loadtxt(uncert_file)
    
    #Create blank lists for reaction numbers and uncertainties
    reaction_numbers = []
    reaction_uncertainties = []
    
    #Break each line of uncert into the reaction and associated uncertainty factor and build the uncertainties
    for uncert_line in uncerts:

        number = int(uncert_line[0])
        uncert = uncert_line[1]

        if number > 12: number += 7 #The H + O2 + M reaction has been split; increase reaction numbers afterwards to compensate

        reaction_numbers += [number - 1] #because uncert.txt is one-indexed and Cantera is zero-indexed
        reaction_uncertainties += [uncert]
    
    #Make the reaction uncertainties list into an array
    reaction_uncertainties = np.array(reaction_uncertainties)
    
    #Intitialize the model chemistry so we have access to the Cantera solution
    mumpce_cantera_model.initialize_chemistry()
    
    #Initialize A-factor uncertaintes
    a_factor_uncertainties = np.ones(mumpce_cantera_model.gas.n_reactions,dtype=np.float64) * 2 # uncertainty of 2 is the default
    a_factor_uncertainties[0:31] = np.repeat(1.2,31) #1.2 is the default for the H2 submodel
    a_factor_uncertainties[reaction_numbers] = reaction_uncertainties #Replace defaults with the information from uncertainty_file
    
    #Create the blank parameter_uncertainties array
    parameter_uncertainties = np.zeros(len(mumpce_cantera_model.model_parameter_info))
    
    #Iterate over all available parameters in the model
    for param_number,param_info in enumerate(mumpce_cantera_model.model_parameter_info):
        #Get this parameter's type and corresponding reaction number
        parameter_type = param_info['parameter_type']
        reaction_number = param_info['reaction_number']
        value = param_info['parameter_value']
        
        uncertainty = a_factor_uncertainties[reaction_number]
        if 'A' in parameter_type:
            #For an A-factor, the uncertainty factor is just the number from uncertainty_file
            parameter_uncertainties[param_number] = uncertainty
        if 'E' in parameter_type:
            #For an activation energy, assume that it contributes the same amount to the uncertainty as the A-factor at 1000 K
            #This number is arbitrary
            #value = np.abs(mumpce_cantera_model.get_parameter(param_number))
            parameter_uncertainties[param_number] = min((value + 1000 * ct.gas_constant * np.log(uncertainty))/value,1.2)
    return parameter_uncertainties