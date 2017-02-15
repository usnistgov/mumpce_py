.. mumpce documentation master file, created by
   sphinx-quickstart on Tue Nov  8 13:53:27 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. currentmodule:: mumpce_new

Welcome to MUM-PCE's documentation!
===================================
   
Welcome to the home page for the Method of Uncertainty Minimization using Polynomial Chaos Expansions. This software is a Python package that implements the methodology presented in [C&F paper, PECS paper, JPCA paper]. The software does the following things:

  * Compiles a database of experimental measurements
  * Constrains a physical model against the measurements in the database (optimization)
  * Determines the uncertainty in the physical model parameters based on the uncertainty in the measurements (uncertainty analysis)
  * Identifies measurements that are inconsistent with the constrained model (outlier detection)
  * Identifies measurements that do not strongly constrain the constrained model (experimental design)

This implementation cannot be used out of the box. Instead, it is necessary for the user to create an interface to the user's own code, which will be specific to that application. Two examples of how to do this are provided. One is a toy model which demonstrates how an interface might be written; it is intended to be as complete as possible while also being simple. The other example is an interface to the reaction kinetics program Cantera; this sort of interface probably represents the worst use case possible, with multiple heterogeneous measurements and a highly complex interface to a detailed model.

The package is implemented in a way to be as general as possible, which means that efficiency is often sacrificed in order to implement this generality. Expert users may be able to modify the code in such a way as to make it more efficient for their particular application. No support is provided for this adventure, but please let me know if you are successful.

Contents:
+++++++++

.. toctree::
   :maxdepth: 2
   :includehidden:

   Project
   measurement
   toy
   cantera


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

