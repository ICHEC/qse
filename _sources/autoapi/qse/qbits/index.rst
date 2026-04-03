qse.qbits
=========

.. py:module:: qse.qbits

.. autoapi-nested-parse::

   This module defines the central object in the QSE package: the Qbits object.



Classes
-------

.. autoapisummary::

   qse.qbits.Qbits


Functions
---------

.. autoapisummary::

   qse.qbits.default


Module Contents
---------------

.. py:class:: Qbits(positions=None, labels=None, states=None, cell=None, pbc=None, calculator=None)

   The Qbits object can represent an isolated molecule, or a
   periodically repeated structure.  It has a unit cell and
   there may be periodic boundary conditions along any of the three
   unit cell axes. Information about the qbits (qubit state and
   position) is stored in ndarrays.

   :Parameters: * **positions** (*list of xyz-positions*) -- Qubit positions.  Anything that can be converted to an
                  ndarray of shape (n, 3) will do: [(x1,y1,z1), (x2,y2,z2),
                  ...].
                * **labels** (*list of str*) -- A list of strings corresponding to a label for each qubit.
                * **states** (*list of 2-length arrays.*) -- State of each qubit.
                * **cell** (*qse.Cell | np.ndarray*) -- Unit cell vectors.
                  Either pass a matrix where each row corresponds to a lattice vector
                  or a qse.Cell object.
                  Default value is all zeros.
                * **pbc** (*one or three bool*) -- Periodic boundary conditions flags.  Examples: True,
                  False, 0, 1, (1, 1, 0), (True, False, False).  Default
                  value: False.
                * **calculator** (*calculator object*) -- Used to attach a calculator for doing computation.

   .. rubric:: Examples

   Empty Qbits object:

   .. jupyter-execute::

       import qse
       qs = qse.Qbits()
       print(qs)

   These are equivalent:

   .. jupyter-execute::

       from qse import Qbit, Qbits
       import numpy as np

       # below two methods are equivalent
       a = Qbits(
           positions=np.array([(0, 0, 0), (0, 0, 2)]),
           labels=['qb1', 'qb2'])

       a = Qbits.from_qbit_list(
           [Qbit('qb1', position=(0, 0, 0)), Qbit('qb2', position=(0, 0, 2))])
       print(a)

       xd = np.array([[0, 0, 0], [0.5, 0.5, 0.5]])
       qdim = Qbits(xd)
       qdim.cell = np.eye(3)
       qdim.pbc = True
       qlat = qdim.repeat([3,3,3])
       print(qlat)

   The qdim will have shape = (1,1,1) and qlat will have shape = (3, 3, 3)

   .. note::
       In order to do computation, a calculator object has to attached
       to the qbits object.


   .. py:property:: calc

      Calculator object.


   .. py:property:: pbc

      Reference to pbc-flags for in-place manipulations.


   .. py:method:: set_pbc(pbc)

      Set periodic boundary condition flags.



   .. py:method:: get_pbc()

      Get periodic boundary condition flags.



   .. py:method:: new_array(name, a, dtype=None, shape=None)

      Add new array.

      If *shape* is not *None*, the shape of *a* will be checked.



   .. py:method:: get_array(name, copy=True)

      Get an array.

      Returns a copy unless the optional argument copy is false.



   .. py:method:: set_array(name, a, dtype=None, shape=None)

      Update array.

      If *shape* is not *None*, the shape of *a* will be checked.
      If *a* is *None*, then the array is deleted.



   .. py:method:: copy()

      Return a copy.



   .. py:method:: todict()

      For basic JSON (non-database) support.



   .. py:method:: fromdict(dct)
      :classmethod:


      Rebuild qbits object from dictionary representation (todict).



   .. py:method:: __repr__()

      We use pymatgen type of style to print Qbits class



   .. py:method:: extend(other)

      Extend qbits object by appending qbits from *other*.



   .. py:method:: append(qbit)

      Append qbit to end.



   .. py:method:: __getitem__(indices)

      Return a subset of the qbits.

      :Parameters: **indices** (*int | list | slice*) -- The indices to be returned.

      :returns: *Qbit | Qbits.* -- If indices is a scalar a Qbit object is returned. If indices
                is a list or a slice, a Qbits object with the same cell, pbc, and
                other associated info as the original Qbits object is returned.



   .. py:method:: __delitem__(indices)

      Delete a subset of the qbits.

      :Parameters: **indices** (*int | list*) -- The indices to be deleted.



   .. py:method:: pop(i=-1)

      Remove and return qbit at index *i* (default last).



   .. py:method:: __imul__(m)

      In-place repeat of qbits.



   .. py:method:: draw(radius=None, show_labels=False, colouring=None, units=None, equal_aspect=True)

      Visualize the positions of a set of qubits.

      :Parameters: * **radius** (*float | str*) -- A cutoff radius for visualizing bonds.
                     Pass 'nearest' to set the radius to the smallest
                     distance between the passed qubits.
                     If no value is passed the bonds will not be visualized.
                   * **show_labels** (*bool*) -- Whether to show the labels of the qubits.
                     Defaults to False.
                   * **colouring** (*str | list*) -- A set of integers used to assign different colors to each Qubit.
                     This can be used to view different magnetic orderings.
                     Must have the same length as the number of Qubits.
                   * **units** (*str, optional*) -- The units of distance.
                   * **equal_aspect** (*bool, optional*) -- Whether to have the same scaling for the axes.
                     Defaults to True.

      .. seealso:: :obj:`qse.draw`



   .. py:method:: repeat(rep)

      Create new repeated qbits object.

      The *rep* argument should be a sequence of three positive
      integers like *(2,3,1)* or a single integer (*r*) equivalent
      to *(r,r,r)*.



   .. py:method:: translate(displacement)

      Translate qbit positions.

      :Parameters: **displacement** (*float | np.ndarray*) -- The displacement argument can be a float an xyz vector or an
                   nx3 array (where n is the number of qbits).



   .. py:method:: get_centroid()

      Get the centroid of the positions.

      :returns: *np.ndarray* -- The centroid of the positions.

      .. rubric:: Notes

      For a set of :math:`k` positions
      :math:`\textbf{x}_1, \textbf{x}_2, ..., \textbf{x}_k`
      the centroid is given by

      .. math::

          \frac{\textbf{x}_1 + \textbf{x}_2 + ... + \textbf{x}_k}{k}.



   .. py:method:: set_centroid(centroid)

      Set the centroid of the positions.

      :Parameters: **centroid** (*float | np.ndarray*) -- The new centroid. Can be a float or a xyz vector

      .. rubric:: Notes

      For a set of :math:`k` positions
      :math:`\textbf{x}_1, \textbf{x}_2, ..., \textbf{x}_k`
      the centroid is given by

      .. math::

          \frac{\textbf{x}_1 + \textbf{x}_2 + ... + \textbf{x}_k}{k}.



   .. py:method:: rotate(a, v='z', center=(0, 0, 0), rotate_cell=False)

      Rotate qbits based on a vector and an angle, or two vectors.

      :Parameters: * **a** -- Angle that the qbits is rotated (anticlockwise) around the vector 'v'.
                     'a' can also be a vector and then 'a' is rotated
                     into 'v'. If 'a' is an angle it must be in degrees.
                   * **v** -- Vector to rotate the qbits around. Vectors can be given as
                     strings: 'x', '-x', 'y', ... .
                     Defaults to 'z'.
                   * **center** -- The center is kept fixed under the rotation. Use 'COP' to
                     fix the center of positions or 'COU' to fix the center of
                     cell. Defaults to = (0, 0, 0).
                   * **rotate_cell = False** -- If true the cell is also rotated.

      .. rubric:: Examples

      The following all rotate 90 degrees anticlockwise around the z-axis,
      so that the x-axis is rotated into the y-axis:

      >>> qbits.rotate(90)
      >>> qbits.rotate(90, 'z')
      >>> qbits.rotate(90, (0, 0, 1))
      >>> qbits.rotate(-90, '-z')
      >>> qbits.rotate('x', 'y')
      >>> qbits.rotate((1, 0, 0), (0, 1, 0))

      .. rubric:: Notes

      If 'a' is an angle, :math:`\theta`, and if :math:`\textbf{v}` is the vector
      then we define

      .. math::
          R = \cos(\theta)I + \sin(\theta)[\textbf{v}]_\times
          + (1-\cos(\theta))\textbf{v}\textbf{v}^T

      where :math:`[\textbf{v}]_\times \textbf{x} = \textbf{v} \times \textbf{x}`.
      If
      :math:`\textbf{r}` is a coordinate vector
      and :math:`\textbf{c}` is the center, this transforms
      the coordinate vector to

      .. math::

          \textbf{r} \rightarrow R(\textbf{r}-\textbf{c}) + \textbf{c}.



   .. py:method:: euler_rotate(phi=0.0, theta=0.0, psi=0.0, center=(0, 0, 0))

      Rotate qbits via Euler angles (in degrees).

      :Parameters: * **phi** (*float*) -- The 1st rotation angle around the z axis.
                   * **theta** (*float*) -- Rotation around the x axis.
                   * **psi** (*float*) -- 2nd rotation around the z axis.
                   * **center** -- The point to rotate about. A sequence of length 3 with the
                     coordinates, or 'COM' to select the center of mass, 'COP' to
                     select center of positions or 'COU' to select center of cell.

      .. rubric:: Notes

      Let

      .. math::

          R =
          \begin{pmatrix}
          \cos(\psi ) & \sin(\psi ) & 0 \\
          -\sin(\psi ) & \cos(\psi ) & 0\\
          0 & 0 & 1\\
          \end{pmatrix}
          \begin{pmatrix}
          1 & 0 & 0\\
          0 & \cos(\theta ) & \sin(\theta ) \\
          0 & -\sin(\theta ) & \cos(\theta ) \\
          \end{pmatrix}
          \begin{pmatrix}
          \cos(\phi ) & \sin(\phi ) & 0 \\
          -\sin(\phi ) & \cos(\phi ) & 0\\
          0 & 0 & 1\\
          \end{pmatrix}

      then if :math:`\textbf{r}` is a coordinate vector
      and :math:`\textbf{c}` is the center, this transforms
      the coordinate vector to

      .. math::

          \textbf{r} \rightarrow R(\textbf{r}-\textbf{c}) + \textbf{c}.



   .. py:method:: get_angle(i: int, j: int, k: int)

      Get the angle in degress formed by three qubits.

      :Parameters: * **i** (*int*) -- The index of the first qubit.
                   * **j** (*int*) -- The index of the second qubit.
                   * **k** (*int*) -- The index of the third qubit.

      :returns: *float* -- The angle between the qubits.

      .. rubric:: Notes

      Let x1, x2, x3 be the vectors describing the positions of the three
      qubits. Then we calcule the angle between x1-x2 and x3-x2.



   .. py:method:: get_angles(indices)

      Get the angle in degress formed by three qubits for multiple groupings.

      :Parameters: **indices** (*list | np.ndarray*) -- The indices of the groupings of qubits.
                   Must be of shape (n, 3), where n is the number of groupings.

      :returns: *np.ndarray* -- The angles between the qubits.

      .. rubric:: Notes

      Let x1, x2, x3 be the vectors describing the positions of the three
      qubits. Then we calcule the angle between x1-x2 and x3-x2 for all the
      different groupings.



   .. py:method:: set_angle(a1, a2=None, a3=None, angle=None, mask=None, indices=None, add=False)

      Set angle (in degrees) formed by three qbits.

      Sets the angle between vectors *a2*->*a1* and *a2*->*a3*.

      If *add* is `True`, the angle will be changed by the value given.

      Same usage as in :meth:`ase.Qbits.set_dihedral`.
      If *mask* and *indices*
      are given, *indices* overwrites *mask*. If *mask* and *indices*
      are not set, only *a3* is moved.



   .. py:method:: rattle(stdev=0.001, seed=None, rng=None)

      Randomly displace qbits.

      This method adds random displacements to the qbit positions.
      The random numbers are
      drawn from a normal distribution of standard deviation stdev.

      For a parallel calculation, it is important to use the same
      seed on all processors!



   .. py:method:: get_distance(i, j)

      Return the distance between two qbits.

      :Parameters: * **i** (*int*) -- The index of the first qubit.
                   * **j** (*int*) -- The index of the second qubit.

      :returns: *float* -- The distance between the qubits.



   .. py:method:: get_distances(i, indices)

      Return distances of the ith qubit with a list of qubits.

      :Parameters: * **i** (*int*) -- The index of the ith qubit.
                   * **indices** (*list[int]*) -- The indices of other qubits.

      :returns: *np.ndarray* -- An array containing the distances.



   .. py:method:: get_all_distances()

      Return the distances of all of the qubits with all of the other qubits.

      :returns: *np.ndarray* -- An array of shape (nqbits, nqbits) containing the distances.



   .. py:method:: set_distance(i, j, distance, fix=0.5, mask=None, indices=None, add=False, factor=False)

      Set the distance between qubits i and j.

      :Parameters: * **i** (*int*) -- The index of the ith qubit.
                   * **j** (*int*) -- The index of the jth qubit.
                   * **distance** (*float*) -- The new distance to be set.
                   * **fix** (*float*) -- By default, the center of the two qbits will be fixed.  Use
                     fix=0 to fix the first qbit, fix=1 to fix the second
                     qbit and fix=0.5 (default) to fix the center of the bond.
                   * **mask** (*np.ndarray | list*) -- If mask or indices are set (mask overwrites indices),
                     only the qbits defined there are moved. It is assumed that the
                     qbits in mask/indices move together with the jth qubit.
                     If fix=1, only the ith qubit will therefore be moved.
                   * **indices** (*np.ndarray | list*) -- If mask or indices are set (mask overwrites indices),
                     only the qbits defined there are moved. It is assumed that the
                     qbits in mask/indices move together with the jth qubit.
                     If fix=1, only the ith qubit will therefore be moved.
                   * **add** -- When add is true, the distance is changed by the value given.
                   * **factor** -- When factor is true, the value given is a factor scaling the distance.



   .. py:method:: wrap(**wrap_kw)

      Wrap positions to unit cell.

      Parameters:

      wrap_kw: (keyword=value) pairs
          optional keywords `pbc`, `center`, `pretty_translation`, `eps`,
          see :func:`ase.geometry.wrap_positions`



   .. py:method:: __eq__(other)

      Check for identity of two qbits objects.

      Identity means: same positions, states, unit cell and
      periodic boundary conditions.



   .. py:method:: __ne__(other)

      Check if two qbits objects are not equal.

      Any differences in positions, states, unit cell or
      periodic boundary condtions make qbits objects not equal.



   .. py:method:: write(filename, format=None, **kwargs)

      Write qbits object to a file.

      see ase.io.write for formats.
      kwargs are passed to ase.io.write.



   .. py:method:: compute_interaction_hamiltonian(distance_func, interaction, tol=1e-08)

      Compute the interaction Hamiltonian for a system of qubits.

      This function constructs an Operators object based on the distances between
      the qubits.

      :Parameters: * **distance_func** (*callable*) -- A function that takes a distance (float) and returns the interaction
                     coefficient (float).
                   * **interaction** (*str | list[str]*) -- The type of interaction (e.g., "X", "Y", "Z") for the Hamiltonian terms,
                     or can be a list of strings.
                   * **tol** (*float, optional*) -- Tolerance threshold for including interaction terms. Terms with absolute
                     coefficients less than `tol` are discarded. Default is 1e-8.

      :returns: *Operators* -- The interaction operators.

      .. rubric:: Examples

      To create a ZZ Hamiltonian for only nearest neighbour qubits

      >>> spacing = 1.0
      >>> qbits = qse.lattices.chain(spacing, 4)
      >>> coupling = -2.
      >>> qbits.compute_interaction_hamiltonian(
      ...     lambda x: coupling*np.isclose(x, spacing), "Z"
      ... )
      ... Number of qubits: 4
      ... Number of terms: 3
      ...
      ... -2.00 Z0 Z1
      ... -2.00 Z1 Z2
      ... -2.00 Z2 Z3

      To create an XY Hamiltonian based on distance

      >>> spacing = 1.0
      >>> qbits = qse.lattices.chain(spacing, 2)
      >>> coupling = 1.
      >>> hamiltonian = qbits.compute_interaction_hamiltonian(
      ...     lambda x: coupling / x**3, ["X", "Y"]
      ... )
      >>> hamiltonian += qbits.compute_interaction_hamiltonian(
      ...     lambda x: coupling / x**3, ["Y", "X"]
      ... )
      ... Number of qubits: 2
      ... Number of terms: 2
      ...
      ... 1.00 X0 Y1
      ... 1.00 Y0 X1



   .. py:method:: add_dim()

      Adds a spatial dimension to the positions. E.e. to go from 1D to 2D systems
      or 2D to 3D systems.



   .. py:method:: remove_dim(dim)

      Removes a spatial dimension to the positions. E.e. to go from 2D to 1D systems
      or 3D to 2D systems.

      :Parameters: **dim** (*str*) -- The dimension to be removed.
                   Must be one of 'x', 'y' or 'z'.



   .. py:method:: get_scaled_positions()

      Get the positions in units of the unit cell.

      :returns: *np.ndarray* -- The positions in units of the unit cell.



.. py:function:: default(data, dflt)

   Helper function for setting default values.


