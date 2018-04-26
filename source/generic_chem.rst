Generic chemistry interface for Cantera
***************************************

The generic chemistry interface for Cantera defines a number of common methods that are used for all Cantera-based simulations. These methods create and modify Cantera phase objects for sensitivity analysis. 

Initialization function summary
===============================

.. automodule:: cantera_utils
   :members: measurement_initialize,measurement_initialize_pd

.. currentmodule:: cantera_chemistry_model

Generic chemistry model method summary
======================================

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


Generic Cantera chemistry model class
=====================================

.. autoclass:: CanteraChemistryModel
   
   .. automethod:: CanteraChemistryModel.sensitivity
   .. automethod:: CanteraChemistryModel.get_parameter
   .. automethod:: CanteraChemistryModel.perturb_parameter
   .. automethod:: CanteraChemistryModel.reset_model
   .. automethod:: CanteraChemistryModel.get_model_parameter_info
   .. automethod:: CanteraChemistryModel.prepare_chemistry
   .. automethod:: CanteraChemistryModel.initialize_chemistry
   .. automethod:: CanteraChemistryModel.blank_chemistry
   