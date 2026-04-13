qse.qubo
========

.. py:module:: qse.qubo

.. autoapi-nested-parse::

   QUBO
   ----

   This module helps frame the QUBO (Quadratic Unconstrained
   Binary Optimization) problem. The Qubo problem is finding
   minima of the following function

   .. math::

       F(x_1, x_2, ..., x_N) = \sum_{i, j} Q_{ij} x_i x_j


   where :math:`x_i` are binary variables, :math:`i` runs over
   the number of variables :math:`N` and the matrix :math:`Q_{ij}`
   defines the QUBO problem.

   Several interesting usecases map to QUBO, such as

   - Ising model
   - Set partitioning problem
   - Portfolio optimization
   - Traffic optimization



Classes
-------

.. autoapisummary::

   qse.qubo.Qubo


Module Contents
---------------

.. py:class:: Qubo(N, N2=None, C6=1.0, Q=None)

   Qubo Class

   :Parameters: * **N** (*int*) -- Problem size, for 1D lattice system,
                  extent in one dimension. For 2D system
                  this takes extent in first dimension
                * **N2** (*int, optional*) -- If provides, the system is 2D, and N2
                  takes the extent in second dimension, by default None
                * **C6** (*float, optional*) -- Coefficient of interaction, by default 1.0
                * **Q** (*np.ndarray, optional*) -- The matrix for the QUBO problem, by default None


   .. py:method:: cost_function(x)

      Cost function to minimize for mapping positions
      for a given Q

      :Parameters: **x** (*np.ndarray*) -- Flattened (original shape :math:`N \times 3` as coordinates)
                   array of length :math:`3N`



