#import mumpce_py as mumpce
import mumpce

#from cantera_chemistry_model import cantera_chemistry_model
from flame_speed import FlameSpeed
#from shock_tube import shock_tube
import shock_tube_utils as stu
import reactions as rxns
#from shock_tube_utils import shock_tube_delay,shock_tube_concentration,shock_tube_ratio

import pandas as pd
import cantera as ct
import numpy as np
import math

def rxn_initialize(name=None,
                   T=None,
                   Patm=None,
                   fuels=None,
                   chemistry_model=None,
                   reaction=None,
                   reaction_denominator=None,
                   measurement_type=None,
                   value=None,
                   uncertainty=None,
                   comment=None,**kwargs
                  ):
    """Initialize a reaction rate constant measurement"""
    
    model_args = [T,Patm,fuels,chemistry_model]
    
    if measurement_type.lower() == 'single':
        model_keys = dict(reaction_number=reaction,**kwargs)
        model = rxns.ReactionRateAtCondition
    if measurement_type.lower() == 'ratio':
        model_keys = dict(reaction_numerator=reaction,reaction_denominator=reaction_denominator,**kwargs)
        model = rxns.ReactionRateRatioAtCondition
    if measurement_type.lower() == 'aratio':
        model_keys = dict(reaction_numerator=reaction,reaction_denominator=reaction_denominator,**kwargs)
        model = rxns.ReactionARatio
    if measurement_type.lower() == 'eratio':
        model_keys = dict(reaction_numerator=reaction,reaction_denominator=reaction_denominator,**kwargs)
        model = rxns.ReactionEDiff
    
    mdl = model(*model_args,**model_keys)
    meas = rxns.RxnMeasurement(name=name,model=mdl,value=value,uncertainty=uncertainty,
                              active_parameters=None,parameter_uncertainties=None,comment=comment,
                              response_type='log'
                             )
    
    return meas
        
    

def ign_initialize(name=None,
                   T=None,
                   Patm=None,
                   fuels=None,
                   initial_timestep=1.0e-5,
                   integration_time=500,
                   critical_value=None,
                   critical_species=None,
                   critical_type=None,
                   chemistry_model=None,
                   critical_denominator=None,
                   critical_rise=None,
                   value=None,
                   uncertainty=None,
                   comment=None,**kwargs
                  ):
    """Initialize a shock tube experiment
    """
    #initial_timestep = 1.0e-4
    fcn = None
    
    kwargs = dict(loglevel=None,**kwargs)
    args = [T,Patm,fuels,ct.IdealGasReactor,chemistry_model]#,fcn)
    
    #if name.startswith('pro') or critical_type == 'tau':
    if critical_type == 'tau':
        #Species concentration at specific time
        integration_time = integration_time/1e6 #Convert from microseconds to seconds        
        model = stu.ShockTubeConcentration
        kwargs = dict(crit_ID=critical_species,integration_time=integration_time,**kwargs)
        
    #elif name.startswith('pul') or critical_type == 'ratio':
    elif critical_type == 'ratio':
        #Concentration ratio at specific time
        integration_time = integration_time/1e6 #Convert from microseconds to seconds
        model = stu.ShockTubeRatio
        kwargs = dict(crit_numerator=critical_species,
                      crit_denom=critical_denominator,
                      integration_time=integration_time,
                      **kwargs)
    else:
        #This is a problem that is finding the delay of something
        model = stu.ShockTubeDelay
        if critical_type == 'crit':
            #Critical species maximum production
            fcn = stu.critical_species_production
            if critical_species == 'PRES':
                #Maximum pressure rise
                fcn = stu.pressure_rise
        if critical_type == 'pres':
            #Maximum pressure rise
            fcn = stu.pressure_rise
        if critical_type == 'conc':
            fcn = stu.target_concentration
        args += [fcn]
        if value: #An experimental value has been supplied. Calculate initial timestep based on this initial value
            base = 10 #hard-coded value, this seems to work well
            delay = math.exp(value)
            #print(delay)
            log_optimal_timestep = math.floor(math.log(delay,base)) - 1
            initial_timestep = (base ** log_optimal_timestep)/1.0e6
        kwargs = dict(crit_ID=critical_species,initial_timestep=initial_timestep,
                      critical_rise=critical_rise,critical_value=critical_value**kwargs)
    
    mdl = model(*args,**kwargs)
    meas = mumpce.Measurement(name=name,model=mdl,value=value,uncertainty=uncertainty,
                              active_parameters=None,parameter_uncertainties=None,comment=comment,
                              response_type='log'
                             )
    
    return meas

def fls_initialize(name=None,
                   T=None,
                   Patm=None,
                   fuels=None,
                   chemistry_model=None,
                   value=None,
                   uncertainty=None,
                   comment=None,**kwargs
                  ):
    mdl = FlameSpeed(T,Patm,fuels,chemistry_model,domain_length=2,initial_points=20,loglevel=0,name=name,**kwargs)
    meas = mumpce.Measurement(name=name,model=mdl,value=value,uncertainty=uncertainty,
                              active_parameters=None,parameter_uncertainties=None,comment=comment)
    return meas
def measurement_initialize(filename,chemistry_model):
    """Read a text file to initialize a batch of measurements
    """
    #Define some parameters so things are easy to find
    num_fuels = 3
    fuel_position = 2
    ox_position = 8
    dil_name_position = 9
    temp_position = 10
    pres_position = temp_position + 1
    crit_species_position = pres_position + 2
    crit_type_position = pres_position + 3
    crit_value_position = pres_position + 4
    crit_denom_position = pres_position + 5
    exp_value_position = pres_position + 7
    exp_uncert_position = pres_position + 6
    
    measurement_list = []
    
    #Open the input file for reading
    with open(filename,'r') as input_file:
        for input_line in input_file:
            
            #Open the single-line input file for reading
            #input_file = open(filename,'r') #Open the master input file
            #input_file.seek(0) #Make sure we're at the beginning of the file
                
            #Read and split the input line
            #input_line = input_file.readline()
            parsed_input = input_line.split()
            
            exp_type = parsed_input[0]
            
            name = parsed_input[0] + '_' + parsed_input[1]
            
            temperature = float(parsed_input[temp_position])
            pressure = float(parsed_input[pres_position])
            critical_species = parsed_input[crit_species_position]
            critical_type = parsed_input[crit_type_position]
            if parsed_input[crit_value_position] == '-':
                critical_value = None
            else:
                critical_value = float(parsed_input[crit_value_position])/1.0e6
            if parsed_input[crit_denom_position] == '-':
                critical_denominator = None
            else:
                critical_denominator = parsed_input[crit_denom_position]
            
            #Build the composition string
            fuel_string = ''
            fuel_total = float(0.0)
            for i in range(num_fuels):
                pos = fuel_position + 2*i
                
                fuel_name = parsed_input[pos]
                fuel_conc = parsed_input[pos+1]
                if float(fuel_conc) > 0:
                    if i > 0:
                        fuel_string = fuel_string + ','
                    fuel_string = fuel_string + fuel_name + ":" + fuel_conc
                    fuel_total = fuel_total + float(fuel_conc)
            #Oxidizer
            ox_conc = parsed_input[ox_position]
            fuel_string = fuel_string + ',O2:' + ox_conc
    
            #Bath gas
            dil_name = parsed_input[dil_name_position]
            if dil_name == 'Air':
                dil_conc = 3.76 * float(ox_conc)
                fuel_string = fuel_string + ',' + 'N2' +':' + str(dil_conc)
            else:
                dil_conc = 1.0 - float(fuel_total) - float(ox_conc)
                fuel_string = fuel_string + ',' + dil_name +':' + str(dil_conc)
            
            #Build the experimental dictionary
            
            print (name)
            
            if exp_type == 'fls':
                meas = fls_initialize(name=name,
                                      T=temperature,
                                      Patm=pressure,
                                      fuels=fuel_string,
                                      chemistry_model=chemistry_model
                                     )
            else:#if exp_type == 'ign':
                meas = ign_initialize(name=name,
                                      T=temperature,
                                      Patm=pressure,
                                      fuels=fuel_string,
                                      critical_species=critical_species,
                                      critical_type=critical_type,
                                      chemistry_model=chemistry_model,
                                      critical_value=critical_value,
                                      critical_denominator=critical_denominator
                                     )
            
                

            
            meas.value = float(parsed_input[exp_value_position])
            meas.uncertainty = float(parsed_input[exp_uncert_position])
            measurement_list += [meas]
    
    return measurement_list#temperature_list,pressure_list,fuel_string_list,critical_species_list

def measurement_initialize_xl(filename,chemistry_model=None,**kwargs):
    """Read a database file in Excel into a Pandas dataframe, then process the dataframe into a batch of measurements
    
    :param filename: The file that contains the experimental database.
    :key chemistry_model: The Cantera chemistry model. It must be a chemistry model that can be used to make a Cantera phase object
    :type filename: str
    """
    
    measurement_list = measurement_initialize_pd(filename,chemistry_model=chemistry_model,**kwargs)
    
    return measurement_list

def measurement_initialize_pd(source,chemistry_model=None,**kwargs):
    """Read a database file in Excel into a Pandas dataframe, then process the dataframe into a batch of measurements
    
    :param source: The file that contains the experimental database.
    :key chemistry_model: The Cantera chemistry model. It must be a chemistry model that can be used to make a Cantera phase object
    :type filename: str
    """
    #Create the blank measurement list
    measurement_list = []
    
    #Read the Excel data for this project
    if type(source) == type(''):
        #This is a string representing a filename, so read it from Excel into Pandas
        df = pd.read_excel(source)
    else:
        #This is a Pandas datafram containing the data
        df = source
    
    #df = pd.read_excel(filename)
    df_columns = df.columns.values
    
    #Find the columns that define the fuel species (this will be in order)
    fuel_names = [s for s in df_columns if 'Fuel' in s]
    fuel_moles = [s for s in df_columns if 'X_fuel' in s]
    
   
    
    #critical_time = None
    critical_denominator = None
    critical_value = None
    critical_rise = None
    reaction_denominator = None
    reaction_numerator = None
    time = None
    value = None
    uncertainty = None
    chemistry = None
    ox = None
    dil = None
    commen = None
    #Check to see if critical_value exists (not all experiments will have one)
    if 'Crit_val' in df_columns:
        critical_value=True
    
    #Find the temperature, pressure, simulation name, critical species, and critical value keywords
    for s in df_columns:
        if s.startswith('Temp'): temp_keyw = s #Temperature
        if s.startswith('Pres'): pres_keyw = s #Pressure
        if s.startswith('Sim'): sim_keyw = s #Simulation name
        if s.startswith('Time'): #Integration time
            time_keyw = s 
            time = True
        if s.startswith('Ox'): #Oxygen specified
            ox_keyw = s
            ox = True
        if s.startswith('Dil'): #Diluent specified
            dil_keyw = s
            dil = True
        if s.startswith('Reac'): #This is a Reaction keyword
            if 'denom' in s:
                rxn_denom_keyw = s
                reaction_denominator = True
            else:
                rxn_keyw = s
                reaction_numerator = True
        if s.startswith('Crit'): #This is a Critical keyword
            if 'spec' in s: #Critical species
                crit_spec_keyw = s
            if 'val' in s:
                crit_val_keyw = s
                critical_value = True
            if 'denom' in s:
                crit_denom_keyw = s
                critical_denominator = True
            if 'rise' in s:
                crit_rise_keyw = s
                critical_rise = True
        if s.startswith('Exp'): #Experimental measurement data
            if 'val' in s:
                exp_val_keyw = s
                value = True
            if 'unc' in s:
                exp_unc_keyw = s
                uncertainty = True
        if s.startswith('Com'): #This project includes comments
            comment_keyw = s
            commen = True
            
    
    #If there is a measurement value, must have uncertainty
    if value is True and uncertainty is False:
        raise ValueError
    
    #Check to see if the database specifies models per experiment (not all projects will need it)
    if 'Model' in df_columns:
        chemistry= True
    
    #Group the experiments by unique ID so that we can iterate over them
    gb = df.groupby('ID')
    for nm,this_experiment in gb:    
        #Get the name of this experiment
        name = this_experiment.Type.values[0] + '_' + this_experiment.ID.values[0]
        #Get the fuels and fuel mole fractions for this experiment
        this_fuels = this_experiment[fuel_names].values[0]
        this_moles = this_experiment[fuel_moles].values[0]
        
        #Parse the fuels information and build the Cantera composition string
        fuel_string = ''
        fuel_total = float(0.0)
        for fuel,moles in zip(this_fuels,this_moles):
            #Don't add a fuel if the mole fraction is 0
            if float(moles) > 0:
                #Needs a comma on the second and subsequent fuel
                if len(fuel_string) > 0:
                    fuel_string = fuel_string + ','
                fuel_string = fuel_string + fuel + ':' + str(moles)
                fuel_total = fuel_total + moles
    
        #Oxidizer
        if ox:
            if not np.isnan(this_experiment.Ox.values[0]): #If the ox value is blank, this experiment doesn't use it
                fuel_string = fuel_string + ',O2:' + str(this_experiment.Ox.values[0])
        
        #Bath gas
        if dil:
            # Air is a special case - Assume that there are 3.76 moles of N2 per mole of oxygen
            if this_experiment[dil_keyw].values[0] == 'Air':
                if ox is None: raise ValueError # Raise an error if Ox is not specified
                dil_concentration = 3.76 * this_experiment[ox_keyw].values[0]
                fuel_string = fuel_string + ',' + 'N2' +':' + str(dil_concentration)
            # Special case for air with argon instead of nitrogen
            elif this_experiment[dil_keyw].values[0] == 'ArAir':
                if ox is None: raise ValueError # Raise an error if Ox is not specified
                dil_concentration = 3.76 * this_experiment[ox_keyw].values[0]
                fuel_string = fuel_string + ',' + 'AR' +':' + str(dil_concentration)
            # If there is any other bath gas specified, it is assumed that the mole fractions sum to 1 and the bath gas is "whatever is left"
            else:
                try:
                    if np.isnan(this_experiment[dil_keyw].values[0]): #If the dil value is blank, don't do anything
                        pass
                except TypeError: #np.isnan() will return an error if you try to give it a string, which means that Dil is not blank
                    dil_concentration = 1.0 - float(fuel_total) - this_experiment[ox_keyw].values[0]
                    fuel_string = fuel_string + ',' + this_experiment[dil_keyw].values[0] +':' + str(dil_concentration)
    
        #chem_ext = this_experiment.Model.values[0]
        if chemistry:
            this_chem = this_experiment.Model.values[0]
            try:
                if np.isnan(this_chem): #Check to see if a chemistry model is specified 
                    chem = chemistry_model
            except TypeError: #np.isnan() will return an error if you try to give it a string, which means that Model is not blank
                chem = this_chem
        else:
            chem = chemistry_model            
        
        #print (chem)
        
        val = None
        unc = None
        if value:
            val = this_experiment[exp_val_keyw].values[0]
            unc = this_experiment[exp_unc_keyw].values[0]
        
        comment = None
        if commen:
            this_comment = this_experiment[comment_keyw].values[0]
            try:
                if np.isnan(this_comment):
                    comment = ''
            except TypeError:
                comment = this_comment
        else:
            comment = ''
            
        
        if this_experiment.Type.values[0] == 'fls':
            meas = fls_initialize(name=name,value=val,uncertainty=unc,
                                  T=this_experiment[temp_keyw].values[0],
                                  Patm=this_experiment[pres_keyw].values[0],
                                  fuels=fuel_string,
                                  chemistry_model=chem,
                                  comment=comment,**kwargs
                                 )
        elif this_experiment.Type.values[0] == 'rxn':
            denom = None
            if reaction_denominator:
                denom = this_experiment[rxn_denom_keyw].values[0]
                if not np.isnan(denom):
                    denom = int(denom)
            meas = rxn_initialize(name=name,value=val,uncertainty=unc,
                                  T=this_experiment[temp_keyw].values[0],
                                  Patm=this_experiment[pres_keyw].values[0],
                                  fuels=fuel_string,
                                  chemistry_model=chem,
                                  measurement_type=this_experiment[sim_keyw].values[0],
                                  reaction=int(this_experiment[rxn_keyw].values[0]),
                                  reaction_denominator=denom,
                                  comment=comment,**kwargs
                                 )
        else:
            cv = None
            rise = None
            denom = None
            tim = None
            if time:
                tim = this_experiment[time_keyw].values[0]
                if tim < 1.0:
                    print('Specified integration time is very small ({} microseconds). ')
                    raise ValueError
            if critical_value:
                cv = this_experiment[crit_val_keyw].values[0]
            if critical_denominator:
                denom = this_experiment[crit_denom_keyw].values[0]
            if critical_rise:
                rise = this_experiment[crit_rise_keyw].values[0]
            meas = ign_initialize(name=name,value=val,uncertainty=unc,
                                  T=this_experiment[temp_keyw].values[0],
                                  Patm=this_experiment[pres_keyw].values[0],
                                  fuels=fuel_string,
                                  critical_species=this_experiment[crit_spec_keyw].values[0],
                                  critical_type=this_experiment[sim_keyw].values[0],
                                  chemistry_model=chem,
                                  critical_value=cv,
                                  integration_time=tim,
                                  critical_denominator=denom,
                                  critical_rise=rise,
                                  comment=comment,**kwargs
                                 )
        
        measurement_list += [meas]
    return measurement_list