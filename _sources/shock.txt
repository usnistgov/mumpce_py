Shock tube models
*****************

.. contents::

This class of model uses a Cantera :class:`cantera.Reactor`, by default a constant volume reactor, as a surrogate for a shock tube. The :py:class:`.StateDefinition` is used as the initial state of the reactor.

Base shock tube class
=====================

.. currentmodule:: shock_tube_base
.. autoclass:: ShockTube
   
   .. automethod:: initialize_reactor
   .. autoinstanceattribute:: reactor_model
   .. automethod:: reset_model
   

Shock tube delay
================

.. currentmodule:: shock_tube_utils
.. autoclass:: ShockTubeDelay
   
   .. automethod:: evaluate
   .. automethod:: optimal_timestep
   .. autoinstanceattribute:: initial_timestep
   
   
   
Generic critical evaluation function
++++++++++++++++++++++++++++++++++++

.. autofunction:: generic_critical_function   

Shock tube concentration
========================

.. autoclass:: ShockTubeConcentration
   :members:

Shock tube ratio
================

.. autoclass:: ShockTubeRatio
   :members:
