#from shock_tube_base import shock_tube
import shock_tube_base as stb
import numpy as np
import cantera as ct
import copy
import math

class ShockTubeDelay(stb.ShockTube):
    """A model for determining the delay time in shock tube ignition
    
    This class provides a model for calculating a shock tube ignition delay time in Cantera. It is based off of the tutorial. 
    
    This class implements the  :func:`evaluate` method required by :func:`mumpce.model`. The :func:`evaluate` method of this class will find the delay time based on one of three possible conditions. The critical_function argument will determine which of these three conditions is being tracked.
    
         * :func:`critical_species_production` - Time at which the critical species production rate is maximized
         * :func:`target_concentration` - Time at which the critical species concentration reaches a certain value
         * :func:`pressure rise` - Time at which the rate of pressure rise is maximized
    
    This class takes all the arguments of the shock_tube class and the following arguments:
    
    :param critical_function: The function that will determine whether the ignition event has occured. Several options are provided
    :keyword crit_ID: The name of the species being tracked. If critical_function is pressure_rise, this does not need to be specified.
    :keyword critical_value: The value of the critical species mole fraction at which integration will stop, if the ignition event is defined as the critical species reaching a certain concentration. Not used otherwise.
    :keyword critical_rise: Whether the delay event is defined by the critical species rising above or falling below critical_value
    :keyword initial_timestep: The initial timestep for integrating the reactor. Default 10 microseconds
    :type critical_function: function
    :type crit_ID: str
    :type critical_value: float
    :type critical_rise: str
    :type initial_timestep: float
    
    """
    def __init__(self,
                 T,Patm,composition,
                 reactor_model,chemistry_model,
                 critical_function,crit_ID=None,critical_value=None,critical_rise=None,
                 initial_timestep=1.0e-5,loglevel=None):
        
        super(ShockTubeDelay,self).__init__(T,Patm,composition,reactor_model,chemistry_model,loglevel)
        
        self.critical_ID = crit_ID
        #self.critical_denominator = crit_denom
        self.critical = critical_function 
        self.critical_value = critical_value
        self.critical_rise = critical_rise
        
        self.initial_timestep = initial_timestep #: The initial timestep used to find the ignition delay. See :py:func:`optimal_timestep`
        
    def __str__(self):
        
        if self.critical_value is None:
            str_args = (self.initial.T,
                        self.initial.P/1.0e3,
                        self.initial.composition,
                        self.critical_ID,
                       )
            modelstr = 'Ignition delay time: {:8.0f} K, {:5.2f} kPa, {}, delay based on maximum d{}/dt'.format(*str_args)
        else:
            str_args = (self.initial.T,
                        self.initial.P/1.0e3,
                        self.initial.composition,
                        self.critical_ID,
                        critical_value*1.0e6,
                       )
            modelstr = 'Ignition delay time: {:8.0f} K, {:5.2f} kPa, {}, delay until mole fraction of {} = {} '.format(*str_args)
        
        #modelstr = 'Ignition delay time: {:8.0d} K, {:5.2d} kPa, {}, delay based on maximum d{}/dt'.format(*str_args)
        
        #modelstr = 'Ignition delay time: ' + str(self.initial.T) + ' K, ' + str(self.initial.P) + ' Pa ' + str(self.initial.composition) + ' ' + self.critical_ID
        return modelstr
        
    def find_delay(self,time_so_far,timestep):
        time = 0.0               
        break_loop = True
        numsteps = 1
        
        critical_last = 0
        
        #Initialize rollback variables
        reactor_contents = copy.deepcopy(self.reactor.thermo.X)
        reactor_temperature = copy.deepcopy(self.reactor.thermo.T)
        reactor_pressure = copy.deepcopy(self.reactor.thermo.P)
        
        contents_interim = copy.deepcopy(self.reactor.thermo.X)
        temperature_interim = copy.deepcopy(self.reactor.thermo.T)
        pressure_interim = copy.deepcopy(self.reactor.thermo.P)
        
        #Integration loop
        while break_loop:
            #keep_going,crit = self.critical(self.reactor,self.critical_ID,critical_last)
            keep_going,crit = self.critical(self,critical_last)
            
            if time_so_far < 1.0e-5:
                keep_going = True
            
            critical_last = crit
            #Print the reactor state
            if not(self.loglevel is None):
                print ('%10.4e %10.4e %10.1f %10.1f %10.8e %10.3e' % (time_so_far,time, 
                                                                      self.reactor.thermo.T, self.reactor.thermo.P, 
                                                                      crit,timestep
                                                                     )
                      )
            
            new_time = time + timestep
            numsteps += 1
            
            if keep_going:
                #The stopping criterion has not been reached
                time = new_time
                time_so_far += timestep
                
                #Save reactor state for restoration purposes
                reactor_contents = copy.deepcopy(contents_interim)#copy.deepcopy(reactor.thermo.X)
                reactor_temperature = copy.deepcopy(temperature_interim)#copy.deepcopy(reactor.thermo.T)
                reactor_pressure = copy.deepcopy(pressure_interim)#copy.deepcopy(reactor.thermo.P)
                
                contents_interim = copy.deepcopy(self.reactor.thermo.X)
                temperature_interim = copy.deepcopy(self.reactor.thermo.T)
                pressure_interim = copy.deepcopy(self.reactor.thermo.P)
                
                #Advance the reactor
                self.simulation.advance(time)
            else:
                time = 0
                time_so_far -= 2*timestep
                time_so_far = max(time_so_far,0)
                break_loop = False
                if not(self.loglevel is None):
                    print '--'
        
        return time_so_far,reactor_temperature,reactor_pressure,reactor_contents
    
    def evaluate(self):
        """Finds the ignition delay time
        
        Computes the ignition delay time with an iterative procedure. First, the delay time is found to within a precision of self.initial_timestep. Then, the reactor moves back two timesteps, halves the size of the timestep, and repeats until the ignition delay time is found to within one part in :math:`10^-5` of the initial timestep. If the default timestep is :math:`10^-5` s, then the final precision is :math:`10^-10` s
        
        If you know that the initial timestep was not set wisely, you can automatically compute one using self.optimal_timestep().
        
        :returns: Ignition delay time in microseconds
        :rtype: float
        """
        #Initialize the chemistry if it has been blanked
        if self.gas is None:
            self.initialize_chemistry()
        self.gas.TPX = self.initial.T, self.initial.P, self.initial.composition
        #Initialize the reactor
        self.initialize_reactor()
        
        #Determine the timestep to start and the desired precision
        timestep = self.initial_timestep
        precision = timestep * 1.0e-5
        
        #Start at zero time
        time_so_far = 0
        while timestep > precision:
            #Calculate the delay time
            delay,current_T,current_P,current_X = self.find_delay(time_so_far,timestep)
            
            #Set the reactor to the conditions at which the calculation will be restarted and re-initialize the reactor
            self.gas.TPX = current_T,current_P,current_X
            self.initialize_reactor()
            #reac1,shock1,gas = shock_tube_initialize(current_T,current_P,current_X,chemistry_model)
            
            timestep = timestep / 2 # Shrink the timestep
            time_so_far = delay #save the elapsed time 
        
        delay = float(delay / 1.0e-6) # Convert from seconds to microseconds
        
        return delay
    
    def optimal_timestep(self):
        """Compute an optimal initial timestep for this measurement, as the default may be too large or too small. Saves the result in self.inital_timestep"""
        delay = self.evaluate()
        
        log_optimal_timestep = math.floor(math.log(delay,10)) - 1
        
        self.initial_timestep = (10 ** log_optimal_timestep)/1.0e6
        return
class ShockTubeConcentration(stb.ShockTube):
    """A model for determining the concentration of a particular species after a certain integration time.
    
    This class provides a model for calculating the concentration in the shock tube at a specified time.
    
    This class implements the  :func:`evaluate` method required by :func:`mumpce.model`. This class takes all the arguments of the shock_tube class and the following additional arguments:
    
    :param crit_ID: The name of the species being tracked.
    :param integration_time: The time over which the reactor will be integrated. Default 100 microseconds
    :type crit_ID: str
    :type integration_time: float
    
    """
    def __init__(self,
                 T,Patm,composition,
                 reactor_model,chemistry_model,crit_ID,
                 integration_time=1.0e-4,loglevel=None):
        
        super(ShockTubeConcentration,self).__init__(T,Patm,composition,reactor_model,chemistry_model,loglevel)
        
        self.critical_ID = crit_ID
        
        self.integration_time = integration_time #: The time for which the reactor will be integrated
    
    def __str__(self):
        
        str_args = (self.initial.T,                    
                    self.initial.P/1.0e3,
                    self.initial.composition,
                    self.critical_ID,
                    self.integration_time * 1.0e3,
                   )
        modelstr = 'Specified mole frac: {:8.0f} K, {:5.2f} kPa, {}, mole fraction of {} at {} ms'.format(*str_args)
        
        #modelstr = 'Specified mole frac: ' + str(self.initial.T) + ' K, ' + str(self.initial.P) + ' Pa ' + str(self.initial.composition) + ', ' + self.critical_ID + ' at ' + self.integration_time + ' seconds'
        return modelstr
    
    def evaluate(self):
        """Calculates the concentration of the critical species at the specified integration time
        
        :returns: Critical species mole fraction
        :rtype: float
        """
        self.initialize_chemistry()
        self.initialize_reactor()
        
        time = self.integration_time
        
        #num_timesteps = 10
        #times = 
        
        
        #if not(self.loglevel is None):
        #        print ('%10.4e %10.1f %10.1f %10.8e %10.3e' % (time_so_far,time, 
        #                                                              self.reactor.thermo.T, self.reactor.thermo.P, 
        #                                                              crit,timestep
        #                                                             )
        #              )
        
        self.simulation.advance(time)

        crit_X = float(self.reactor.thermo[self.critical_ID].X)
        
        return crit_X[0]
    
class ShockTubeRatio(stb.ShockTube):
    """A model for determining the mole fration ratios of two particular species after a certain integration time.
    
    This class provides a model for calculating the ratio of the concentrations of two particular species in the shock tube at a specified time. The ratio is defined as :math:`y = C_n/C_d` where :math:`C_n` and :math:`C_d` are the concentrations of the numerator and denominator species respectively.
    
    This class implements the  :func:`evaluate` method required by :func:`mumpce.model`. It class takes all the arguments of the shock_tube class and the following additional arguments:
    
    :param crit_numerator: The name of the numerator species 
    :param crit_denom: The name of the denominator species 
    :param integration_time: The time over which the reactor will be integrated. Default 100 microseconds
    :type crit_numerator: str
    :type crit_denom: str
    :type integration_time: float
    
    """
    def __init__(self,
                 T,Patm,composition,
                 reactor_model,chemistry_model,
                 crit_numerator=None,crit_denom=None,integration_time=1.0e-4,loglevel=None):
        
        super(ShockTubeRatio,self).__init__(T,Patm,composition,reactor_model,chemistry_model,loglevel)
        
        self.critical_numerator = crit_numerator
        self.critical_denominator = crit_denom
        
        self.integration_time = integration_time #: The time for which the reactor will be integrated
        
    def __str__(self):
        
        str_args = (self.initial.T,                    
                    self.initial.P/1.0e3,
                    self.initial.composition,
                    self.critical_numerator,
                    self.critical_denominator,
                    self.integration_time * 1.0e3
                   )
        
        modelstr = 'Concentration ratio: {:8.0f} K, {:5.2f} kPa, {}, [{}]/[]{} at {} ms'.format(*str_args)
        
#        modelstr = 'Concentration ratio: ' + str(self.initial.T) + ' K, ' + str(self.initial.P) + ' Pa ' + str(self.initial.composition) + ', [' + self.critical_numerator +']/[' + self.critical_denominator + '] at ' + self.integration_time + ' seconds'
        return modelstr
    def evaluate(self):
        """Compute the concentration of the critical species
        
        :returns: Critical species mole fraction
        :rtype: float
        """
        self.initialize_chemistry()
        self.initialize_reactor()
        
        time = self.integration_time
        
        self.simulation.advance(time)

        crit_numerator = self.reactor.thermo[self.critical_numerator].X
        crit_denominator = self.reactor.thermo[self.critical_denominator].X
        
        ratio = float(crit_numerator / crit_denominator)
        
        return ratio[0]
    
    
def generic_critical_function(measurement,critical_last):
    """A generic function to define whether the ignition delay criterion has been satisfied for use with the :func:`shock_tube_delay` class
    
    This is a function used by :py:class:`ShockTubeDelay` to determine if the ignition delay criterion has been satisfied. If it has been satisified, the integration will stop. :py:class:`ShockTubeDelay` stores this 
    The function is called in the following way: ::
    
         keep_going,crit = self.critical(self,critical_last)
    
    If keep_going is True, :py:class:`ShockTubeDelay` continues integrating. If False, it will stop. 
    
    When crit is calculated, it can be compared with critical_last. critical_last is output from the critical_function at the last timestep.
    
    :param measurement: The shock tube object
    :param critical_last: The crit output from the previous timestep.
    :returns: keep_going: Whether the integration should continue
    :rtype: keep_going bool
    :returns: crit
    :rtype: crit float
    
    """
    
    return keep_going,crit
def critical_species_production(measurement,crit_rop_last):
    
    crit_spec = measurement.critical_ID
    crit_index = measurement.reactor.thermo.kinetics_species_index(crit_spec)
    
    #crit_rop_last = crit_rop
    crit_create = measurement.reactor.thermo.creation_rates[crit_index]
    crit_destro = measurement.reactor.thermo.destruction_rates[crit_index]
    crit_rop = crit_create - crit_destro
    
    crit_rop_max = max(crit_rop_last,crit_rop)
    
    keep_going = crit_rop >= crit_rop_max
    
    #print crit_rop, crit_rop_last
    
    return keep_going, crit_rop

def pressure_rise(measurement,pressure_rise_last):
    #roc = shock_tube.thermo.creation_rates
    #rod = shock_tube.thermo.destruction_rates
    rop = measurement.reactor.thermo.net_production_rates
    
    energies = measurement.reactor.thermo.partial_molar_int_energies
    temp,density = measurement.reactor.thermo.TD
    temp,pressure = measurement.reactor.thermo.TP
    cv_mass = measurement.reactor.thermo.cv_mass
    molwt = measurement.reactor.thermo.mean_molecular_weight
    total_concentration = pressure / (ct.gas_constant * temp)
    dconcentration_dt = rop.sum()
    
    # dT/dt = -1 * sum_k(molar_energy*molar_production_rate) / (density * cv_mass)
    enthalpy_rate = -1 * np.dot(rop, energies.T) / (density * cv_mass)
    #print enthalpy_rate
    
    # dP/dt = gas_constant * (concentration*dT/dt + temp*dConc/dt
    
    pressure_rise = ct.gas_constant * (total_concentration * enthalpy_rate + dconcentration_dt * temp)

    pressure_rise = max(0,pressure_rise)
    
    pressure_rise_max = max(pressure_rise,pressure_rise_last)
    
    keep_going = pressure_rise >= pressure_rise_max
    
    return keep_going, pressure_rise

def target_concentration(measurement,concentration_last):
    crit_concentration = measurement.reactor.thermo[crit_spec].X
    crit_val = measurement.critical_value
    
    if measurement.critical_rise == 'rise':
        keep_going = crit_concentration < crit_val
    else:
        keep_going = crit_concentration > crit_val
        
    return keep_going,crit_concentration