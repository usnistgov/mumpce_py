Projects and Solutions
**********************

MUM-PCE is logically divided into projects. A project is represented by a :py:class:`.Project` object which contains the following items:

   * A list of :py:class:`.Measurement` objects which define the experimental database that will constrain the model parameters
   * A separate list of :py:class:`.Measurement` objects which define the application space of the model
   * A *parametric model* which defines the uncertain parameters that will be constrained by the measurements
   * A list of *physics models* implemented by :py:class:`.Model` objects which define how to use the parametric model to perform a calculation for each :py:class:`.Measurement`
   * An :py:func:`.initialize_function` that will read outside data 
   * A :py:class:`.Solution` object which contains the constrained parameters and a structure that defines their uncertainty

The intent of this implementation is that the user only ever interacts with the :py:class:`.Project` object. The user has to create the :py:class:`Model` objects and also a suitable initialization function. Once this is done, however, :py:class:`.Project` has the necessary tools to evaluate the measurements, optimize the model, and perform the outlier and experimental design analyses.

Creating Projects
=================

A blank :py:class:`.Project` can be instantiated with no arguments::
   
   my_project = mumpce.Project()

This creates a :py:class:`.Project` with no measurement list, no intialization function, and no physics model of any kind. Such a project would not be very useful, and so normally the user would define a model object and some measurement initialization function::
   
   my_project = mumpce.Project(initialize_function=my_initialize,
                               model=my_model,
                               parameter_uncertainties=my_uncertainty)

Initializing Measurements within the Project
++++++++++++++++++++++++++++++++++++++++++++

A :py:class:`.Project` that has been created with an initialization function, as above, allows the initialization function to be accessed directly from the :py:class:`.Project`::
   
   my_project.measurement_initialize('file_with_my_measurements.ext')

This call has the initialization function read the file 'file_with_my_measurements.ext' and then populate the measurement list.

Creating a Project with existing Measurements
+++++++++++++++++++++++++++++++++++++++++++++

Alternatively, the user can initialize the measurement list outside of the Project and then include that list as an argument when creating the Project::
   
   my_measurements = my_initialize('file_with_my_measurements.ext',my_model)
   my_applications = my_initialize('file_with_my_applications.ext',my_model)
   my_project = mumpce.Project(measurement_list=my_measurements,
                               application_list=my_applications,
                               parameter_uncertainties=my_uncertainty)

Solving Projects to create Solutions
====================================

A Project that has a measurement list and parameter uncertainties is ready to be solved. The workflow is:

   * Perform sensitivity analysis on the measurements::
      
      my_project.find_sensitivity()
      
   * Find parameters with sensitivity greater than a cutoff value and define those as active::
      
      my_project.find_active_parameters(sens_cutoff)
      
   * Pass the active parameters to the measurements in the measurement list::
      
      my_project.set_active_parameters()
      
   * Calculate the response function for each measurement in the measurement list, used for the optimization::
      
      my_project.make_response()
      
   * Solve for the constrained model based on the response functions::
      
      my_project.run_optimization()
      
   * Identify and remove the inconsistent measurements::
      
      my_project.remove_inconsistent_measurements()
      
   * Identify and remove the measurements that do not constrain the model::
      
      my_project.remove_low_information_measurements()

This process will produce a :py:class:`.ResponseSurface` object within each :py:class:`.Measurement` and a :py:class:`Solution` object associated with the Project as a whole.

:py:class:`.Project` method summary
===================================

.. currentmodule:: Project

.. autosummary::
   Project
   Project.measurement_initialize
   Project.application_initialize
   Project.find_sensitivity
   Project.find_active_parameters
   Project.set_active_parameters
   Project.make_response
   Project.run_optimization
   Project.validate_solution
   Project.remove_inconsistent_measurements
   Project.calculate_entropy
   Project.remove_low_information_measurements
   Project.plot_pdfs
   
Project class
=============
.. currentmodule:: Project
.. autoclass:: Project

   .. autoinstanceattribute:: measurement_list
   .. autoinstanceattribute:: initialize_function
   .. automethod:: Project.measurement_initialize
   .. autoinstanceattribute:: application_list
   .. autoinstanceattribute:: app_initialize_function
   .. automethod:: Project.application_initialize
   .. automethod:: Project.find_active_parameters
   .. automethod:: Project.run_optimization
   .. automethod:: Project.validate_solution
   .. automethod:: Project.calculate_entropy
   .. automethod:: Project.plot_pdfs
   
   

Solution class
==============

.. currentmodule:: solution
.. autoclass:: Solution
   
   .. autoinstanceattribute:: Solution.x
   .. autoinstanceattribute:: Solution.cov
   .. autoinstanceattribute:: Solution.alpha
   