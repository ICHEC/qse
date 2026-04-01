qse.signal
==========

.. py:module:: qse.signal

.. autoapi-nested-parse::

   Definition of the Signal class.



Classes
-------

.. autoapisummary::

   qse.signal.Signal


Module Contents
---------------

.. py:class:: Signal(values, duration=None)

   The Signal class represents a 1D signal with values and duration.

   :Parameters: * **values** (*list | np.ndarray*) -- The values of the signal.
                * **duration** (*int*) -- Duration of the signal.
                  Defaults to the length of the passed values.

   .. rubric:: Examples

   One can create a signal by passing an array of values:

   >>> qse.Signal([1, 2, 3], 100)
   ... Signal(duration=100, values=[1. 2. 3.])

   Arithmetic operations with scalars is supported.
   Adding or multiplying a scalar to a Signal returns
   a new Signal with modified values and the same duration.
   For example:

   >>> signal = qse.Signal([1, 1])
   >>> signal * 3 + 0.5
   ... Signal(duration=2, values=[3.5 3.5])

   Two Signals can be added together which
   concatenates their values and sums their durations.
   So if w1, w2 are instantiation of Signal, then
   w = w1 + w2 gives signal with concatenated values, i.e.,
   w.values = [w1.values, w2.values], and added duration
   w.duration = w1.duration + w2.duration.
   For example:

   >>> signal_1 = qse.Signal([1, 1], 5)
   >>> signal_2 = qse.Signal([2, 2], 2)
   >>> signal_1 + signal_2
   ... Signal(duration=7, values=[1. 1. 2. 2.])

   .. rubric:: Notes

   Currently, the object gets created for multi-dim arrays as well.
   However, it should be used for 1D only, we haven't made it useful
   or consistent for multi-dim usage.


   .. py:method:: __iter__()

      Iterate over the signal values.

      :returns: *iterator* -- An iterator over the values.



   .. py:method:: __getitem__(i)

      Access individual signal value(s) by index or slice.

      :Parameters: **i** (*int or slice*) -- Index or slice to access.

      :returns: *scalar or ndarray* -- The corresponding value(s).



   .. py:method:: __eq__(other) -> bool

      Check for equality with another signal.

      :Parameters: **other** (*Signal*) -- Signal to compare with.

      :returns: *bool* -- True if equal in both values and duration.



   .. py:method:: __add__(other)

      Add another signal or scalar to this signal.

      :Parameters: **other** (*Signal or float or int*) -- Signal to concatenate or scalar to add elementwise.

      :returns: *Signal* -- The resulting signal.

      :raises TypeError: If the operand type is unsupported.



   .. py:method:: __radd__(other)

      Right addition operator.

      :returns: *Signal* -- The resulting signal.



   .. py:method:: __iadd__(other)

      In-place addition with another signal or scalar.

      :Parameters: **other** (*Signal or float or int*) -- Signal to concatenate or scalar to add elementwise.

      :returns: *Signal* -- The updated signal.

      :raises TypeError: If the operand type is unsupported.



   .. py:method:: __mul__(other)

      Multiply signal values by a scalar.

      :Parameters: **other** (*float or int*) -- Scalar multiplier.

      :returns: *Signal* -- The resulting signal.

      :raises TypeError: If the operand type is unsupported.



   .. py:method:: __rmul__(other)

      Right multiplication by a scalar.

      :returns: *Signal* -- The resulting signal.



   .. py:method:: __imul__(other)

      In-place multiplication by a scalar.

      :Parameters: **other** (*float or int*) -- Scalar multiplier.

      :returns: *Signal* -- The updated signal.

      :raises TypeError: If the operand type is unsupported.



   .. py:method:: __repr__() -> str

      Return a string representation of the signal.

      :returns: *str* -- String representation.



