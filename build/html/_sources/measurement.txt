Measurements and Models
***********************

Measurements
============

A :py:class:`.Measurement` represents a physical measurement. The :py:class:`.Measurement` object contains the following items:
   
   * A :py:class:`.Model` object that tells the :py:class:`.Measurement` how to simulate the actual physical measurement. This returns a single number. 
   * The *parametric model* that contains the uncertain parameters.
   * A :py:class:`.ResponseSurface` object that contains information how the py:class:`.Model` varies with respect to the parameters in the parametric model.

The logic behind MUM-PCE is that a measurement is a single number. Calling the :py:func:`evaluate` method of the :py:class:`.Measurement` will call :py:class:`.Model` one time and return that number. This can be inefficient, especially if the :py:class:`.Model` objects actually contain many useful pieces of information.

Initializing measurements
+++++++++++++++++++++++++

The initialization function must be written by the user and will be specific to the user's application. The toy model example :py:func:`.toy_initialize` and the Cantera interface :py:func:`.measurement_initialize_pd`  provide examples for how a user might do this. Both examples use Pandas to read an Excel spreadsheet containing the experimental database. Each line of the spreadsheet must contain enough information to completely describe each simulation that is being performed. The initialization function reads this spreadsheet, creates a :py:class:`.Model` object to simulate the measurement, and then a :py:class:`.Measurement` object containing that :py:class:`.Model`.

Models and Response Surfaces
============================

Defining models
+++++++++++++++

The :py:class:`.Model` object is the interface between MUM-PCE and the user's code The :py:class:`.Model` class provided with MUM-PCE is a template that the user can follow to interface their own code with MUM-PCE. It is defined as a with a set of abstract methods, and the intent is that the user's own model class will be a subclass. The abstract methods must be overwritten by a subclass, or else there will be an error. 

In order to be a valid model, the user's model class must provide an :py:func:`evaluate` method that will return the model values :math:`y_i` and a :py:func:`sensitivity` method that will return the sensitivity coefficients :math:`S_ij = d\ln y_i / d\ln x_j`. In addition, the class must provide methods that can retrive or perturb the parameters within the parametric model and also a list of dictionary-like containers that explain what the parameters are. More information is available in the documentation of :py:class:`.Model`.

Response Surfaces
+++++++++++++++++

The :py:class:`.ResponseSurface` object is the structure that the Project actually interacts with when it is calculating the constrained model. Its interface mimics that of :py:class:`.Model`, insofar as it has an :py:func:`evaluate` and :py:func:`sensitivity` method, which returns more or less the same information. This object will be created by 

Function summary
================

:py:class:`.Measurement` method summary
+++++++++++++++++++++++++++++++++++++++

.. currentmodule:: measurement
.. autosummary::
   Measurement
   Measurement.__str__
   Measurement.evaluate
   Measurement.evaluate_sensitivity
   Measurement.make_response
   Measurement.evaluate_response
   Measurement.sensitivity_response
   Measurement.evaluate_uncertainty
   Measurement.save
   Measurement.load

:py:class:`.Model` method summary
+++++++++++++++++++++++++++++++++

.. currentmodule:: model
.. autosummary:: 
   Model
   Model.__str__
   Model.evaluate
   Model.sensitivity
   Model.get_parameter
   Model.perturb_parameter
   Model.get_model_parameter_info

:py:class:`.ResponseSurface` method summary
++++++++++++++++++++++++++++++++++++++++++++

.. currentmodule:: response_surface
.. autosummary::
   ResponseSurface
   ResponseSurface.evaluate
   ResponseSurface.sensitivity
   


Measurement class
=================

.. currentmodule:: measurement
.. autoclass:: Measurement
   
   .. automethod:: Measurement.evaluate
   .. automethod:: Measurement.evaluate_sensitivity
   .. automethod:: Measurement.make_response
   .. automethod:: Measurement.evaluate_response
   .. automethod:: Measurement.sensitivity_response
   .. automethod:: Measurement.evaluate_uncertainty
   .. automethod:: Measurement.save
   .. automethod:: Measurement.load



Generic model class
===================

.. currentmodule:: model
.. autoclass:: Model
   
   .. automethod:: Model.evaluate
   .. automethod:: Model.sensitivity
   .. automethod:: Model.get_parameter
   .. automethod:: Model.perturb_parameter
   .. automethod:: Model.get_model_parameter_info


Response surface class
======================

.. currentmodule:: response_surface

.. autoclass:: ResponseSurface
   
   .. autoinstanceattribute:: z
   .. autoinstanceattribute:: a
   .. autoinstanceattribute:: b

   .. automethod:: ResponseSurface.evaluate
   .. automethod:: ResponseSurface.sensitivity 