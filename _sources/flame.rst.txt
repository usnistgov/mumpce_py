Laminar flame speeds
********************

This model uses a Cantera :class:`cantera.FreeFlame` to determine the propagation speed of a freely-propagating laminar flame. 

Laminar flame speed model
=========================

.. currentmodule:: flame_speed

.. autoclass:: FlameSpeed
   
   .. automethod:: evaluate
   .. automethod:: initialize_reactor
   .. automethod:: save_restart
   .. automethod:: load_restart
   .. automethod:: ignore_restart