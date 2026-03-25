QSE
===

The Quantum Simulation Environment (QSE) package is intended to provide a flexible,
modular way to frame a quantum simulation problem involving position dependend quantum degrees of freedom.
Primarily we look at a collection of qubits at given set of coordinates.
These can be arbitrary coordinates, or defined on a lattice.

QSE's design is adapted from Atomic Simulation Environment (ASE) to suit the needs
for an abstract representation for

#. `defining quantum computing systems`
#. `computing operations/simulations`

in a vendor agnostic way.
ASE's modular nature, and extensability make it very useful for a similar quantum computing application.

.. mermaid::
   graph TB

   A[QSE]
   Qs[Qbits]
   calc[Calculator]

   A ----> Qs & calc


.. autoclasstree:: qse.Qbits


.. note::

   This project is under active development.


.. code-block:: python

   import qse
   qsqr = qse.lattices.square(
      lattice_spacing=2.0,
      repeats_x=4, repeats_y=4)
   qsqr.draw(radius=5.0)



Contributing
------------
See the `contributing page <https://github.com/ICHEC/qse/blob/main/CONTRIBUTIONS.md>`_.


.. toctree::
   :hidden:
   :maxdepth: 1

   Home <self>
   INSTALLATION

.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: Tutorials

   tutorials/creating_and_manipulating_qbits
   tutorials/generating_lattices
   tutorials/pulser-calc-example
   tutorials/ssh_model


.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: API Reference

   api 
