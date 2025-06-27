qse.qbits
=========

.. py:module:: qse.qbits

.. autoapi-nested-parse::

   Definition of the Qbits class.

   This module defines the central object in the QSE package: the Qbits object.



Classes
-------

.. autoapisummary::

   qse.qbits.Qbits


Functions
---------

.. autoapisummary::

   qse.qbits.string2vector
   qse.qbits.default


Module Contents
---------------

.. py:class:: Qbits(labels=None, states=None, positions=None, scaled_positions=None, cell=None, pbc=None, celldisp=None, constraint=None, calculator=None, info=None)

   The Qbits object can represent an isolated molecule, or a
   periodically repeated structure.  It has a unit cell and
   there may be periodic boundary conditions along any of the three
   unit cell axes. Information about the qbits (qubit state and
   position) is stored in ndarrays.

   :Parameters: * **labels** (*list of str*) -- A list of strings corresponding to a label for each qubit.
                * **states** (*list of 2-length arrays.*) -- State of each qubit.
                * **positions** (*list of xyz-positions*) -- Qubit positions.  Anything that can be converted to an
                  ndarray of shape (n, 3) will do: [(x1,y1,z1), (x2,y2,z2),
                  ...].
                * **scaled_positions** (*list of scaled-positions*) -- Like positions, but given in units of the unit cell.
                  Can not be set at the same time as positions.
                * **cell** (*3x3 matrix or length 3 or 6 vector*) -- Unit cell vectors.  Can also be given as just three
                  numbers for orthorhombic cells, or 6 numbers, where
                  first three are lengths of unit cell vectors, and the
                  other three are angles between them (in degrees), in following order:
                  [len(a), len(b), len(c), angle(b,c), angle(a,c), angle(a,b)].
                  First vector will lie in x-direction, second in xy-plane,
                  and the third one in z-positive subspace.
                  Default value: [0, 0, 0].
                * **celldisp** (*Vector*) -- Unit cell displacement vector. To visualize a displaced cell
                  around the center of mass of a Systems of qbits. Default value
                  = (0,0,0)
                * **pbc** (*one or three bool*) -- Periodic boundary conditions flags.  Examples: True,
                  False, 0, 1, (1, 1, 0), (True, False, False).  Default
                  value: False.
                * **constraint** (*constraint object(s)*) -- Used for applying one or more constraints during structure
                  optimization.
                * **calculator** (*calculator object*) -- Used to attach a calculator for doing computation.
                * **info** (*dict of key-value pairs*) -- Dictionary of key-value pairs with additional information
                  about the system.  The following keys may be used by ase:

                    - spacegroup: Spacegroup instance
                    - unit_cell: 'conventional' | 'primitive' | int | 3 ints
                    - adsorbate_info: Information about special adsorption sites

                  Items in the info attribute survives copy and slicing and can
                  be stored in and retrieved from trajectory files given that the
                  key is a string, the value is JSON-compatible and, if the value is a
                  user-defined object, its base class is importable.  One should
                  not make any assumptions about the existence of keys.

   .. rubric:: Examples

   Empty Qbits object:

   >>> qs = qse.Qbits()

   These are equivalent:

   >>> a = qse.Qbits(
   ...     labels=['qb1', 'qb2'],
   ...     positions=np.array([(0, 0, 0), (0, 0, 2)])
   ... )
   >>> a = qse.Qbits.from_qbit_list(
   ...     [Qbit('qb1', position=(0, 0, 0)), Qbit('qb2', position=(0, 0, 2))]
   ... )

   >>> xd = np.array(
   ...    [[0, 0, 0],
   ...     [0.5, 0.5, 0.5]])
   >>> qdim = qse.Qbits(positions=xd)
   >>> qdim.cell = [1,1,1]
   >>> qdim.pbc = True
   >>> qlat = qdim.repeat([3,3,3])

   The qdim will have shape = (2,1,1) and qlat will have shape = (6, 3, 3)

   .. rubric:: Notes

   In order to do computation, a calculator object has to attached
   to the qbits object.


   .. py:property:: calc

      Calculator object.


   .. py:property:: shape

      The shape of the qbits


   .. py:method:: set_constraint(constraint=None)

      Apply one or more constrains.

      The *constraint* argument must be one constraint object or a
      list of constraint objects.



   .. py:method:: set_cell(cell, scale_qbits=False, apply_constraint=True)

      Set unit cell vectors.

      Parameters:

      cell: 3x3 matrix or length 3 or 6 vector
          Unit cell.  A 3x3 matrix (the three unit cell vectors) or
          just three numbers for an orthorhombic cell. Another option is
          6 numbers, which describes unit cell with lengths of unit cell
          vectors and with angles between them (in degrees), in following
          order: [len(a), len(b), len(c), angle(b,c), angle(a,c),
          angle(a,b)].  First vector will lie in x-direction, second in
          xy-plane, and the third one in z-positive subspace.
      scale_qbits: bool
          Fix qbit positions or move qbits with the unit cell?
          Default behavior is to *not* move the qbits (scale_qbits=False).
      apply_constraint: bool
          Whether to apply constraints to the given cell.

      Examples:

      Two equivalent ways to define an orthorhombic cell:

      >>> qbits = Qbits('He')
      >>> a, b, c = 7, 7.5, 8
      >>> qbits.set_cell([a, b, c])
      >>> qbits.set_cell([(a, 0, 0), (0, b, 0), (0, 0, c)])

      FCC unit cell:

      >>> qbits.set_cell([(0, b, b), (b, 0, b), (b, b, 0)])

      Hexagonal unit cell:

      >>> qbits.set_cell([a, a, c, 90, 90, 120])

      Rhombohedral unit cell:

      >>> alpha = 77
      >>> qbits.set_cell([a, a, a, alpha, alpha, alpha])



   .. py:method:: set_celldisp(celldisp)

      Set the unit cell displacement vectors.



   .. py:method:: get_celldisp()

      Get the unit cell displacement vectors.



   .. py:method:: get_cell(complete=False)

      Get the three unit cell vectors as a `class`:ase.cell.Cell` object.

      The Cell object resembles a 3x3 ndarray, and cell[i, j]
      is the jth Cartesian coordinate of the ith cell vector.



   .. py:method:: get_cell_lengths_and_angles()

      Get unit cell parameters. Sequence of 6 numbers.

      First three are unit cell vector lengths and second three
      are angles between them::

          [len(a), len(b), len(c), angle(b,c), angle(a,c), angle(a,b)]

      in degrees.



   .. py:method:: get_reciprocal_cell()

      Get the three reciprocal lattice vectors as a 3x3 ndarray.

      Note that the commonly used factor of 2 pi for Fourier
      transforms is not included here.



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



   .. py:method:: has(name)

      Check for existence of array.

      name must be one of: 'momenta', 'masses', 'initial_magmoms',
      'initial_charges'.



   .. py:method:: set_positions(newpositions, apply_constraint=True)

      Set positions, honoring any constraints. To ignore constraints,
      use *apply_constraint=False*.



   .. py:method:: get_positions(wrap=False, **wrap_kw)

      Get array of positions.

      Parameters:

      wrap: bool
          wrap qbits back to the cell before returning positions
      wrap_kw: (keyword=value) pairs
          optional keywords `pbc`, `center`, `pretty_translation`, `eps`,
          see :func:`ase.geometry.wrap_positions`



   .. py:method:: get_properties(properties)

      This method is experimental; currently for internal use.



   .. py:method:: copy()

      Return a copy.



   .. py:method:: todict()

      For basic JSON (non-database) support.



   .. py:method:: fromdict(dct)
      :classmethod:


      Rebuild qbits object from dictionary representation (todict).



   .. py:method:: extend(other)

      Extend qbits object by appending qbits from *other*.



   .. py:method:: append(qbit)

      Append qbit to end.



   .. py:method:: pop(i=-1)

      Remove and return qbit at index *i* (default last).



   .. py:method:: repeat(rep)

      Create new repeated qbits object.

      The *rep* argument should be a sequence of three positive
      integers like *(2,3,1)* or a single integer (*r*) equivalent
      to *(r,r,r)*.



   .. py:method:: translate(displacement)

      Translate qbit positions.

      :Parameters: **displacement** (*float | np.ndarray*) -- The displacement argument can be a float an xyz vector or an
                   nx3 array (where n is the number of qbits).



   .. py:method:: center_in_unit_cell(vacuum=None, axis=(0, 1, 2), about=None)

      Center qbits in unit cell.

      Centers the qbits in the unit cell, so there is the same
      amount of vacuum on all sides.

      vacuum: float (default: None)
          If specified adjust the amount of vacuum when centering.
          If vacuum=10.0 there will thus be 10 Angstrom of vacuum
          on each side.
      axis: int or sequence of ints
          Axis or axes to act on.  Default: Act on all axes.
      about: float or array (default: None)
          If specified, center the qbits about <about>.
          I.e., about=(0., 0., 0.) (or just "about=0.", interpreted
          identically), to center about the origin.



   .. py:method:: get_centroid(scaled=False)

              Get the centroid of the positions.

              Parameters
              ----------
              scaled : bool
                  If scaled=True the centroid in scaled coordinates is returned.

              Notes
              -----
              For a set of $k$ positions $    extbf{x}_1,     extbf{x}_2, ...,        extbf{x}_k$
              the centroid is given by
              $
      rac{  extbf{x}_1 +    extbf{x}_2 + ... +      extbf{x}_k}{k}.$




   .. py:method:: set_centroid(centroid, scaled=False)

              Set the centroid of the positions.

              Parameters
              ----------
              centroid : float | np.ndarray
                  The new centroid. Can be a float or a xyz vector
              scaled : bool
                  If scaled=True the centroid is expected in scaled coordinates.

              Notes
              -----
              For a set of $k$ positions $    extbf{x}_1,     extbf{x}_2, ...,        extbf{x}_k$
              the centroid is given by
              $
      rac{  extbf{x}_1 +    extbf{x}_2 + ... +      extbf{x}_k}{k}.$




   .. py:method:: rotate(a, v, center=(0, 0, 0), rotate_cell=False)

      Rotate qbits based on a vector and an angle, or two vectors.

      :Parameters: * **a** -- Angle that the qbits is rotated around the vector 'v'. 'a'
                     can also be a vector and then 'a' is rotated
                     into 'v'.
                   * **v** -- Vector to rotate the qbits around. Vectors can be given as
                     strings: 'x', '-x', 'y', ... .
                   * **center** -- The center is kept fixed under the rotation. Use 'COP' to
                     fix the center of positions or 'COU' to fix the center of
                     cell. Defaults to = (0, 0, 0).
                   * **rotate_cell = False** -- If true the cell is also rotated.

      .. rubric:: Examples

      Rotate 90 degrees around the z-axis, so that the x-axis is
      rotated into the y-axis:

      >>> qbits = Qbits()
      >>> qbits.rotate(90, 'z')
      >>> qbits.rotate(90, (0, 0, 1))
      >>> qbits.rotate(-90, '-z')
      >>> qbits.rotate('x', 'y')
      >>> qbits.rotate((1, 0, 0), (0, 1, 0))



   .. py:method:: euler_rotate(phi=0.0, theta=0.0, psi=0.0, center=(0, 0, 0))

      Rotate qbits via Euler angles (in degrees).

      See e.g http://mathworld.wolfram.com/EulerAngles.html for explanation.

      :Parameters: * **phi** (*float*) -- The 1st rotation angle around the z axis.
                   * **theta** (*float*) -- Rotation around the x axis.
                   * **psi** (*float*) -- 2nd rotation around the z axis.
                   * **center** -- The point to rotate about. A sequence of length 3 with the
                     coordinates, or 'COM' to select the center of mass, 'COP' to
                     select center of positions or 'COU' to select center of cell.



   .. py:method:: get_dihedral(a0, a1, a2, a3, mic=False)

      Calculate dihedral angle.

      Calculate dihedral angle (in degrees) between the vectors a0->a1
      and a2->a3.

      Use mic=True to use the Minimum Image Convention and calculate the
      angle across periodic boundaries.



   .. py:method:: get_dihedrals(indices, mic=False)

      Calculate dihedral angles.

      Calculate dihedral angles (in degrees) between the list of vectors
      a0->a1 and a2->a3, where a0, a1, a2 and a3 are in each row of indices.

      Use mic=True to use the Minimum Image Convention and calculate the
      angles across periodic boundaries.



   .. py:method:: set_dihedral(a1, a2, a3, a4, angle, mask=None, indices=None)

      Set the dihedral angle (degrees) between vectors a1->a2 and
      a3->a4 by changing the qbit indexed by a4.

      If mask is not None, all the qbits described in mask
      (read: the entire subgroup) are moved. Alternatively to the mask,
      the indices of the qbits to be rotated can be supplied. If both
      *mask* and *indices* are given, *indices* overwrites *mask*.

      **Important**: If *mask* or *indices* is given and does not contain
      *a4*, *a4* will NOT be moved. In most cases you therefore want
      to include *a4* in *mask*/*indices*.

      Example: the following defines a very crude
      ethane-like molecule and twists one half of it by 30 degrees.

      >>> qbits = Qbits('HHCCHH', [[-1, 1, 0], [-1, -1, 0], [0, 0, 0],
      ...                          [1, 0, 0], [2, 1, 0], [2, -1, 0]])
      >>> qbits.set_dihedral(1, 2, 3, 4, 210, mask=[0, 0, 0, 1, 1, 1])



   .. py:method:: rotate_dihedral(a1, a2, a3, a4, angle=None, mask=None, indices=None)

      Rotate dihedral angle.

      Same usage as in :meth:`ase.Qbits.set_dihedral`: Rotate a group by a
      predefined dihedral angle, starting from its current configuration.



   .. py:method:: get_angle(index_1: int, index_2: int, index_3: int, mic: bool = False)

      Get the angle in degress formed by three qbits.

      :Parameters: * **index_1** (*int*) -- The index of the first qubit.
                   * **index_2** (*int*) -- The index of the second qubit.
                   * **index_3** (*int*) -- The index of the third qubit.
                   * **mic** (*bool*) -- Use mic=True to use the Minimum Image Convention and calculate the
                     angle across periodic boundaries.

      .. rubric:: Notes

      Let x1, x2, x3 be the vectors describing the positions of the three
      qubits. Then we calcule the angle between x1-x2 and x3-x2.



   .. py:method:: get_angles(indices, mic=False)

      Get angle formed by three qbits for multiple groupings.

      Calculate angle in degrees between vectors between qbits a2->a1
      and a2->a3, where a1, a2, and a3 are in each row of indices.

      Use mic=True to use the Minimum Image Convention and calculate
      the angle across periodic boundaries.



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

      This method adds random displacements to the qbit positions,
      taking a possible constraint into account.  The random numbers are
      drawn from a normal distribution of standard deviation stdev.

      For a parallel calculation, it is important to use the same
      seed on all processors!



   .. py:method:: get_distance(a0, a1, mic=False, vector=False)

      Return distance between two qbits.

      Use mic=True to use the Minimum Image Convention.
      vector=True gives the distance vector (from a0 to a1).



   .. py:method:: get_distances(a, indices, mic=False, vector=False)

      Return distances of qbit No.i with a list of qbits.

      Use mic=True to use the Minimum Image Convention.
      vector=True gives the distance vector (from a to self[indices]).



   .. py:method:: get_all_distances(mic=False, vector=False)

      Return distances of all of the qbits with all of the qbits.

      Use mic=True to use the Minimum Image Convention.



   .. py:method:: set_distance(a0, a1, distance, fix=0.5, mic=False, mask=None, indices=None, add=False, factor=False)

      Set the distance between two qbits.

      Set the distance between qbits *a0* and *a1* to *distance*.
      By default, the center of the two qbits will be fixed.  Use
      *fix=0* to fix the first qbit, *fix=1* to fix the second
      qbit and *fix=0.5* (default) to fix the center of the bond.

      If *mask* or *indices* are set (*mask* overwrites *indices*),
      only the qbits defined there are moved
      (see :meth:`ase.Qbits.set_dihedral`).

      When *add* is true, the distance is changed by the value given.
      In combination
      with *factor* True, the value given is a factor scaling the distance.

      It is assumed that the qbits in *mask*/*indices* move together
      with *a1*. If *fix=1*, only *a0* will therefore be moved.



   .. py:method:: get_scaled_positions(wrap=True)

      Get positions relative to unit cell.

      If wrap is True, qbits outside the unit cell will be wrapped into
      the cell in those directions with periodic boundary conditions
      so that the scaled coordinates are between zero and one.

      If any cell vectors are zero, the corresponding coordinates
      are evaluated as if the cell were completed using
      ``cell.complete()``.  This means coordinates will be Cartesian
      as long as the non-zero cell vectors span a Cartesian axis or
      plane.



   .. py:method:: set_scaled_positions(scaled)

      Set positions relative to unit cell.



   .. py:method:: wrap(**wrap_kw)

      Wrap positions to unit cell.

      Parameters:

      wrap_kw: (keyword=value) pairs
          optional keywords `pbc`, `center`, `pretty_translation`, `eps`,
          see :func:`ase.geometry.wrap_positions`



   .. py:method:: get_volume()

      Get volume of unit cell.



   .. py:property:: cell

      The :class:`ase.cell.Cell` for direct manipulation.


   .. py:method:: write(filename, format=None, **kwargs)

      Write qbits object to a file.

      see ase.io.write for formats.
      kwargs are passed to ase.io.write.



.. py:function:: string2vector(v)

   Used in rotate method to rotate qbit location


.. py:function:: default(data, dflt)

   Helper function for setting default values.


