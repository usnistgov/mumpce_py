The Cantera interface to MUM-PCE
********************************

The Cantera model is designed to demonstrate a real-world example of MUM-PCE. Like the toy model, it consists of a set of experimental measurements that will be used to constrain a chemistry model.

.. contents::
   :depth: 3

.. toctree::
   :maxdepth: 2
   
   generic_chem
   shock
   flame
   reactions
   cantera_examples


The initialization function
===========================

The function :py:func:`.measurement_initialize_xl` reads an Excel file to determine the experimental database. The header row gives the name of each column, which will be parsed by :py:func:`.measurement_initialize_xl`. The order of the columns does not matter. Two examples are given below. Other examples are given in the repository.

.. tip:: The columns can be in any order. If a column is not used by any measurement in the list, it can be omitted. If a column is not used by a particular measurement, it can be left blank; the parser will ignore it.

Example data file
+++++++++++++++++

Example 1
---------

.. csv-table:: 
   :widths: auto
   :file: Flame speeds example.csv


Example 2
---------

.. csv-table:: 
   :widths: auto
   :file: Ignition delay example.csv

Columns of the data file
++++++++++++++++++++++++

The data file contains the following options. Columns with names marked with an asterisk are required, and all others are optional. Text enclosed in (parentheses) is optional. For instance, the parser will recognize Temp, Temper, and Temperature as valid keywords for specifying the unburned gas temperature.

   * Type*: The type of measurement, which determines which Cantera model will be used
      * ign: Shock tube, :py:class:`.ShockTube` and its subclasses
      * pro: Alias for Sim=tau, :py:class:`.ShockTubeConcentration`
      * pul: Alias for Sim=ratio, :py:class:`.ShockTubeRatio`
      * fls: Laminar flame, :py:class:`.FlameSpeed`
      * rxn: Reaction rate measurement, :py:class:`.ReactionRateBase` and its subclasses
   * ID*: A free-form text identifier
      .. note:: The 'name' of the experiment will be generated from Type and ID
   * Reactant definitions*:
      * Fuel: A text identifier for a reactant species in the unburned gas mixture. This must be a species in the Cantera chemistry model
      * X_Fuel: The mole fraction of each reactant species in the unburned gas mixture. 
         .. tip:: As seen in the example, any number of Fuel and X_fuel pairs can be used. They can be in any order, but the first Fuel column must correspond to the first X_fuel column. Although the columns are named Fuel1, Fuel2, etc in the example, this is arbitrary and any identifier, or none, can be used.
         
      * Ox: The mole fraction of oxygen.
      * Dil: The diluent species. The behavior depends on what this is.
         * 'Air': 3.76 moles of nitrogen is added per mole of oxygen. Ox must be specified
         * 'ArAir': 3.76 moles of argon is added per mole of oxygen. Ox must be specified
         * A species in the Cantera chemistry model: the total mole fractions of all the fuels and oxidizer add to less than unity, then enough of this species will be added to ensure that the mole fractions sum to unity.
         
      .. note:: At least one substance must be specified somehow by one of these keywords. If Dil='Air' or Dil='ArAir', Ox be specified.
   * Temp(erature)*: The temperature of the unburned gas in Kelvins
   * Pres(sure)*: The pressure of the unburned gas in standard atmospheres
   * Model: The Cantera chemistry model used for this measurement. 
      .. tip:: The value of ``chemistry_model`` used in the function call is used as a default for all measurements. If this column is not provided, ``chemistry_model`` is used for all measurements. If this column is provided, a value needs to be specified only for those experiments that use a different model
      
      .. warning:: Different models for different measurements is only intended to be used if some measurements are conducted under conditions where rate constants may change (e.g. certain kinds of pressure dependencies). All models should have exactly the same reactions and corresponding parameters, although they can take different values.
   * Sim(ulation): For :py:class:`.ShockTube` and :py:class:`.ReactionRateBase` measurements, defines which sub-measurement will be done.
      * crit: :py:class:`.ShockTubeDelay`, ignition delay time measured by maximum critical species production
      * pres: :py:class:`.ShockTubeDelay`, ignition delay time measured by maximum pressure rise
      * conc: :py:class:`.ShockTubeDelay`, delay time until a certain critical species concentration
      * tau: :py:class:`.ShockTubeConcentration`, concentration of the critical species after a certain time
      * ratio: Behavior depends on whether this is a :py:class:`.ShockTube` or :py:class:`.ReactionRateBase`
         * :py:class:`.ShockTubeRatio`, ratio of two concentrations after a certain time
         * :py:class:`.ReactionRateRatioAtCondition`, the ratio of two rate constants 
      * single: :py:class:`.ReactionRateAtCondition`, a single rate constant 
      * ARatio: :py:class:`.ReactionARatio`, the ratio of two reaction rate A factors 
      * ERatio: :py:class:`.ReactionEDiff`, The ratio of two reaction rate activation energies 
   * Time: For :py:class:`.ShockTubeConcentration` and :py:class:`.ShockTubeRatio` experiments, this is the integration time in microseconds.
   * Crit(ical\_)spec(ies): The species whose mole fraction is to be tracked. If the value in this column is PRES, then pressure is tracked in the simulation, otherwise it must be a species in the Cantera chemistry model. For :py:class:`.ShockTubeRatio` measurements, this is the species whose concentration appears in the numerator.
      .. note:: For backwards compatibility with previous database implementations, Simulation='crit' and Critical_species='PRES' is equivalent to Simulation='pres'
   * Crit(ical\_)denom(inator): For :py:class:`.ShockTubeRatio` measurements, this is the species whose concentration appears in the denominator.
   * Crit(ical\_)rise: For :py:class:`.ShockTubeDelay` measurements with Simulation='conc', this determines if the critical species concentration should be rising or falling through the critical value
   * Crit(ical\_)val(ue):  For :py:class:`.ShockTubeDelay` measurements with Simulation='conc', this is the critical species mole fraction at which the delay occurs.
   * Reac(tion): For :py:class:`.ReactionRateAtCondition` measurements, the reaction whose rate constant is measured. For other measurements based on :py:class:`.ReactionRateBase` experiments, the reaction whose rate constant appears in the numerator.
   * Reac(tion\_)denom(inator): For other measurements based on :py:class:`.ReactionRateBase`, the reaction whose rate constant appears in the denominator.
   * Exp(erimental\_)val(ue): The measured value for this measurement.
   * Exp(erimental\_)unc(ertainty): The uncertainty in the experimentally-measured value
      .. warning:: If Exp_val is present, Exp_unc must also be preseint


Model classes in the Cantera interface
======================================

Generic interface
+++++++++++++++++

The generic chemistry interface for Cantera, :py:class:`.CanteraChemistryModel`, defines a number of common methods that are used for all Cantera-based simulations. These methods create and modify Cantera phase objects for sensitivity analysis. It defines all of the methods required by a :py:class:`.Model` except for :func:`evaluate`.

Laminar flame speeds
++++++++++++++++++++

This model uses a Cantera :class:`cantera.FreeFlame` to determine the propagation speed of a freely-propagating laminar flame. 


Shock tubes
+++++++++++

This class of model uses a Cantera :class:`cantera.Reactor`, by default a constant volume reactor, as a surrogate for a shock tube. The :py:class:`.StateDefinition` is used as the initial state of the reactor.


Shock tube ignition delay time
------------------------------

This model finds the ignition delay time :math:`\tau_{\text{ign}}` by finding the maximum pressure rise or species concentration production rate, by solving:

.. math::
   
   \tau_{\text{ign}} = \text{argmax}_\tau \frac{dX}{dt}\bigg|_{t=\tau}
   
where :math:`X(t)` is either the pressure or the mole fraction of the critical species in the reactor. Specified with Type = 'ign' and Sim = 'crit'.


Shock tube species profile concentration
----------------------------------------

This model finds the time :math:`\tau_{\text{ign}}` at which the critical species concentration reaches a critical value by solving:

.. math::

  \tau_{\text{crit}} = \text{argmin}_\tau (X(t=\tau) - X_{\text{crit}})^2

where :math:`X(t)` is the mole fraction of the critical species in the reactor and :math:`X_{\text{crit}}` is the critical value. There is an additional constraint that :math:`dX/dt` must be either positive or negative when this occurs. Specified with Type = 'ign' and Sim = 'conc'.

Shock tube species profile time
-------------------------------

This model finds the critical species mole fraction :math:`X_{\text{crit}}` by solving:

.. math::
   
   X_{\text{crit}} = X(t=\tau_{\text{crit}})

where :math:`X(t)` is the mole fraction of the critical species in the reactor and :math:`\tau_{\text{crit}}` is the integration time. Specified with Type = 'pro'.

Shock tube species concentration ratio
--------------------------------------

This model finds the critical species mole fraction ratio :math:`r_{\text{crit}}` by solving:

.. math::
   
   r_{\text{crit}} = \frac{X_{\text{num}}(t)}{X_{\text{denom}}(t)}\bigg|_{t=\tau_{\text{crit}}}

where :math:`X_{\text{num}}(t)` and :math:`X_{\text{denom}}(t)` are the mole fractions of the numerator and denominator species in the reactor and :math:`tau_{\text{crit}}` is the integration time. Specified with Type = 'pul'.

Reaction rate measurements
++++++++++++++++++++++++++

This class of model will directly query the Cantera model for rate constant values and optimize the rate paramters to match that value. There is no Cantera :class:`cantera.Reactor` model, but there is a Cantera :class:`cantera.Solution` model that represents the chemistry model.

Reaction rate at specified condition
------------------------------------

This model specifies the measurement of a rate constant at a particular condition, :math:`k(T,P,X)`.

Reaction rate at specified condition
------------------------------------

This model specifies the ratio of two rates at a particular condition, :math:`k_{num}(T,P,X)/k_{den}(T,P,X)`

Reaction rate A-factor ratio
----------------------------

This model specifies the ratio of two reaction rate A factors, regardless of condition, :math:`A_{num}/A_{den}`

Reaction rate activation energy difference
------------------------------------------

This model specifies the difference between two reaction rate activation energies, regardless of condition, :math:`E_{num} - E_{den}`

Function summary
================

Generic chemistry model summary
+++++++++++++++++++++++++++++++

.. currentmodule:: cantera_chemistry_model

.. autosummary:: 
   CanteraChemistryModel
   CanteraChemistryModel.sensitivity
   CanteraChemistryModel.get_parameter
   CanteraChemistryModel.perturb_parameter
   CanteraChemistryModel.reset_model
   CanteraChemistryModel.get_model_parameter_info
   CanteraChemistryModel.prepare_chemistry
   CanteraChemistryModel.initialize_chemistry
   CanteraChemistryModel.blank_chemistry

.. currentmodule:: state_definition

.. autosummary::
   StateDefinition
   
Laminar flame speed summary
+++++++++++++++++++++++++++

.. currentmodule:: flame_speed

.. autosummary::
   FlameSpeed
   FlameSpeed.evaluate
   FlameSpeed.initialize_reactor

Shock tube summary
++++++++++++++++++

.. currentmodule:: shock_tube_base

.. autosummary::
   ShockTube
   ShockTube.initialize_reactor
   ShockTube.reset_model

Shock tube delay summary
------------------------

.. currentmodule:: shock_tube_utils

.. autosummary::
   ShockTubeDelay
   ShockTubeDelay.evaluate

Shock tube concentration summary
--------------------------------

.. autosummary::
   ShockTubeConcentration
   ShockTubeConcentration.evaluate

Shock tube concentration ratio summary
--------------------------------------

.. autosummary::
   ShockTubeRatio
   ShockTubeRatio.evaluate

Reaction rate summary
+++++++++++++++++++++

Base reaction rate class
------------------------

.. currentmodule:: reactions

.. autosummary::
   RxnMeasurement
   RxnMeasurement.make_response
   ReactionRateBase
   ReactionRateBase.sensitivity

Reaction rate subclasses
------------------------

.. autosummary::
   ReactionRateAtCondition
   ReactionRateRatioAtCondition
   ReactionARatio
   ReactionEDiff