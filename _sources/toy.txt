A toy example using MUM-PCE
***************************
  
This toy model is designed to demonstrate the entire workflow for MUM-PCE. It contains a set of measurements and applications that have hard-wired dependencies on a set of parameters, which guarantees that their response surfaces will take a particular form. This model contains measurements that are known to be outliers and also measurements that poorly constrain the model, thus providing an example of the functionality of MUM-PCE. This document explains how the model is developed, how the initialization function interacts with the model to create measurements, and how these entities interact with a MUM-PCE Project.

In addition, :doc:`toy_model_example` is provided which shows the complete workflow of the toy project:

.. toctree::
   :hidden:
   
   toy_model_example
   


Defining models
===============

The models used in the toy example are defined by the equation,

.. math::
    y = A \prod_i k_i p_i


This equation is equivalent to :math:`\ln y = a^\text{T}x + z` where :math:`x_i = \ln p_i`, :math:`a_i = \ln k_i`, and :math:`z = \ln A`. The parameters :math:`p_i` are initially all set to 1, so that the nominal response of the model is equal to :math:`A`. Consequently, the response surface will be identical to model.

A :py:class:`.toy_model` object can be created directly::
   
   experiment_number = 0
   my_model = toy_model(experiment_number)

.. tip:: This usage is not the normal way of calling MUM-PCE. One would normally use an initialization function to create a batch of :py:class:`.Measurement` and :py:class:`.Model` objects at once.

Calling the :func:`my_model.evaluate` method of the toy model returns :math:`y`. Calling :func:`my_model.sensitivity` will return :math:`y` and also the sensitivity vector :math:`S_i = d\ln y/d\ln k_i`, which is equal to :math:`a_i` to within numerical error. :func:`my_model.sensitivity` has the sensitivity perturbation as a parameter, and the smaller this value is, the closer :math:`S` will be to :math:`a`. 

.. code:: python

    my_model.evaluate()




.. parsed-literal::

    1.2214027581601699



.. code:: python

    V,S = my_model.sensitivity(0.0000001)


.. parsed-literal::

    Value =  1.22140275816
    Param  Value+       Value-           Sensitivity
         0  1.22140e+00   1.22140e+00   1.0000e+00
         1  1.22140e+00   1.22140e+00   6.0000e-01
         2  1.22140e+00   1.22140e+00   4.0000e-01
         3  1.22140e+00   1.22140e+00   1.0000e-01
         4  1.22140e+00   1.22140e+00   3.0000e-02
         5  1.22140e+00   1.22140e+00   1.0000e-02
         6  1.22140e+00   1.22140e+00   2.0000e-02

.. note:: The sensitivity output shown above is printed to the log file and also to the screen. 

Calling :func:`my_model.get_parameter_info` returns a list of dictionaries that has information about the model parameters. The keys of this dictionary are the number of the parameter and its name; the names are (imaginatively) Parameter 1 through Parameter 7. As an example::

    my_model.get_model_parameter_info()

.. parsed-literal::

    array([{'parameter_number': 0, 'parameter_name': 'Parameter 1'},
           {'parameter_number': 1, 'parameter_name': 'Parameter 2'},
           {'parameter_number': 2, 'parameter_name': 'Parameter 3'},
           {'parameter_number': 3, 'parameter_name': 'Parameter 4'},
           {'parameter_number': 4, 'parameter_name': 'Parameter 5'},
           {'parameter_number': 5, 'parameter_name': 'Parameter 6'},
           {'parameter_number': 6, 'parameter_name': 'Parameter 7'}], dtype=object)


The initialization function
===========================

The function :py:func:`.toy_initialize` gives an example of how to instantiate measurements. This function reads an Excel file to determine the experimental database. The Excel file for the measurement list looks like the following:

.. parsed-literal::
   Name			Number	Value	Uncertainty
   Experiment 1		0	0.4	0.05
   Experiment 2		1	0.25	0.05
   Experiment 3		2	0.2	0.05
   Experiment 4		3	0.7	0.08
   Experiment 5		4	0.4	0.8

The header row gives the names of each column, which will be parsed by :py:func:`.toy_initialize`. Each row contains the simulation metadata for one measurement along with the relevant measurement data.

   * "Name" defines a unique text identifier for each measurement. This identifier is entirely freeform, but because file names will be derived from this name, all characters used in the name must be allowable by the file system. This name also defines whether :py:class:`.toy_model` or :py:class:`.toy_app` is used.
   * "Value" is the experimentally-measured value with which the simulation output is to be compared. 
   * "Uncertainty" is the uncertainty in the experimentally-measured value.
   * "Number" is the only piece of simulation metadata for the toy model. This number is the experiment_number used when creating the :py:class:`.toy_model` or :py:class:`.toy_app` object.

.. note:: It is allowed to create a :py:class:`.Measurement` without a value or uncertainty, but it is not possible to optimize a model if any members of the measurement list do not have a value or an uncertainty.

The initialization function is called with the Excel experimental database as its argument::
   
   my_measurements = mumpce.toy.toy_initialize('mumpce_toy_experiments.xlsx',mumpce.toy.toy_model)
   my_apps = mumpce.toy.toy_initialize('mumpce_toy_apps.xlsx',mumpce.toy.toy_app)
   
This creates one list which has the measurements and another with the targeted applications.

Creating a Project
==================

In its most simple form, a mumpce :py:class:`.Project` can be created by defining a measurement initialization function and some parameter uncertainties. ::
   
   my_project = mumpce.Project(parameter_uncertainties=mumpce.toy.parameter_uncertainties,
                                initialize_function=mumpce.toy.toy_initialize)

This creates an almost-empty Project object that is ready to have its measurements initialized.


Initializing Measurements within the Project
++++++++++++++++++++++++++++++++++++++++++++

If the :py:class:`.Project` has no measurements to begin with, the measurement and application lists can be created by the :func:`measurement_initialize` method ::
   
   my_project.measurement_initialize('mumpce_toy_experiments.xlsx')
   my_project.application_initialize('mumpce_toy_apps.xlsx')

Creating a Project with existing Measurements
+++++++++++++++++++++++++++++++++++++++++++++

If the measurement and application lists have already been created, the :py:class:`.Project` can then be created with the measurement and application lists as arguments ::

   my_project = mumpce.Project(measurement_list=toy_measurements,
                             application_list=toy_apps,
                             model=mumpce.toy.toy_model,
                             parameter_uncertainties=mumpce.toy.parameter_uncertainties,
                             initialize_function=mumpce.toy.toy_initialize)



Solving the Project
===================

The :py:class:`.Project` is now ready to be solved. The workflow described in :doc:`Project` can be followed from this point::
   
   my_project.find_sensitivity()
   my_project.find_active_parameters(sens_cutoff)
   my_project.set_active_parameters()
   my_project.run_optimization()
   my_project.remove_inconsistent_measurements()
   my_project.remove_low_information_measurements()

Function summary
================

Toy model class summary
+++++++++++++++++++++++

.. currentmodule:: toy

.. autosummary::
   toy_initialize
   toy_model
   toy_app

Initialization function
+++++++++++++++++++++++

.. autofunction:: toy_initialize


Toy Model
+++++++++

:func:`toy_model` method summary
--------------------------------

.. autosummary::
   toy_model
   toy_model.__str__
   toy_model.evaluate
   toy_model.sensitivity
   toy_model.get_parameter
   toy_model.perturb_parameter
   toy_model.reset_model
   toy_model.get_model_parameter_info

.. autoclass:: toy_model

   .. automethod::  toy_model.evaluate
   .. automethod::  toy_model.sensitivity
   .. automethod::  toy_model.get_parameter
   .. automethod::  toy_model.perturb_parameter
   .. automethod::  toy_model.reset_model
   .. automethod::  toy_model.get_model_parameter_info

.. autoclass:: toy_app
   :members:
   
