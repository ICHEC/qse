qse.lattices
============

.. py:module:: qse.lattices

.. autoapi-nested-parse::

   A collection of convenience functions for creating common lattices.
   Qbits class is very flexible to generate arbitrary structures. For
   convenience, helper functions to generate several common types of
   lattices is defined here. See the examples below -


   .. jupyter-execute::

       import qse
       qsqr = qse.lattices.square(
           lattice_spacing=1.0,
           repeats_x=4, repeats_y=4)
       qsqr.draw(radius=2)




Functions
---------

.. autoapisummary::

   qse.lattices.chain
   qse.lattices.hexagonal
   qse.lattices.kagome
   qse.lattices.ring
   qse.lattices.square
   qse.lattices.torus
   qse.lattices.triangular


Package Contents
----------------

.. py:function:: chain(lattice_spacing: float, repeats: int) -> qse.Qbits

   Generate a Qbits object in linear chain geometry.

   :Parameters: * **lattice_spacing** (*float*) -- The lattice spacing.
                * **repeats** (*int*) -- The number of repeats.
                  Must be greater than 1.

   :returns: *Qbits* -- The Qbits lattice.


.. py:function:: hexagonal(lattice_spacing: float, repeats_x: int, repeats_y: int) -> qse.Qbits

   Generate a Qbits object in hexagonal lattice geometry.

   :Parameters: * **lattice_spacing** (*float*) -- The lattice spacing.
                * **repeats_x** (*int*) -- The number of repeats in the x direction.
                  Must be greater than 1.
                * **repeats_y** (*int*) -- The number of repeats in the y direction.
                  Must be greater than 1.

   :returns: *Qbits* -- The Qbits lattice.


.. py:function:: kagome(lattice_spacing: float, repeats_x: int, repeats_y: int) -> qse.Qbits

   Generate a Qbits object in kagome lattice geometry.

   :Parameters: * **lattice_spacing** (*float*) -- The lattice spacing.
                * **repeats_x** (*int*) -- The number of repeats in the x direction.
                  Must be greater than 1.
                * **repeats_y** (*int*) -- The number of repeats in the y direction.
                  Must be greater than 1.

   :returns: *Qbits* -- The Qbits lattice.


.. py:function:: ring(spacing: float, nqbits: int) -> qse.Qbits

   Generate a Qbits object in a ring geometry.

   :Parameters: * **radius** (*float*) -- The spacing between the qubits.
                * **nqbits** (*int*) -- Number of qubits on the ring.

   :returns: *Qbits* -- The Qbits object.


.. py:function:: square(lattice_spacing: float, repeats_x: int, repeats_y: int) -> qse.Qbits

   Generate a Qbits object in square lattice geometry.

   :Parameters: * **lattice_spacing** (*float*) -- The lattice spacing.
                * **repeats_x** (*int*) -- The number of repeats in the x direction.
                  Must be greater than 1.
                * **repeats_y** (*int*) -- The number of repeats in the y direction.
                  Must be greater than 1.

   :returns: *Qbits* -- The Qbits lattice.


.. py:function:: torus(n_outer: int, n_inner: int, inner_radius: float, outer_radius: float) -> qse.Qbits

   Generate a Qbits object in a torus geometry.

   :Parameters: * **n_outer** (*int*) -- Number of points in larger (outer) dimension.
                * **n_inner** (*int*) -- Number of points in smaller (inner) dimension.
                * **inner_radius** (*float*) -- The inner radius.
                * **outer_radius** (*float*) -- The outer radius.

   :returns: *Qbits* -- The Qbits object.


.. py:function:: triangular(lattice_spacing: float, repeats_x: int, repeats_y: int) -> qse.Qbits

   Generate a Qbits object in triangular lattice geometry.

   :Parameters: * **lattice_spacing** (*float*) -- The lattice spacing.
                * **repeats_x** (*int*) -- The number of repeats in the x direction.
                  Must be greater than 1.
                * **repeats_y** (*int*) -- The number of repeats in the y direction.
                  Must be greater than 1.

   :returns: *Qbits* -- The Qbits lattice.


