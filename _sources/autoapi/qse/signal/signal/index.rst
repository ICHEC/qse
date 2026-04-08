qse.signal.signal
=================

.. py:module:: qse.signal.signal

.. autoapi-nested-parse::

   Definition of the Signal class.



Classes
-------

.. autoapisummary::

   qse.signal.signal.Signal


Module Contents
---------------

.. py:class:: Signal(values, duration=None)

   The Signal class represents a 1D signal with values and duration.
   The duration must be divisble by the number of values and each value
   is assume to last for an equal duration.

   :Parameters: * **values** (*list | np.ndarray*) -- The values of the signal.
                * **duration** (*int*) -- Duration of the signal.
                  Defaults to the length of the passed values.

   .. rubric:: Examples

   To create a constant signal, passing a single value:

   .. jupyter-execute::

       import qse
       s = qse.Signal([1], 10)
       print(s)

   To create an arbitrary signal, pass an array
   whose length is equal to the duration:

   .. jupyter-execute::

       import qse
       import numpy as np

       ss = qse.Signal(np.linspace(0, 1, 5), 5)
       print(ss)

   Arithmetic operations with scalars is supported.
   Adding or multiplying a scalar to a Signal returns
   a new Signal with modified values and the same duration.
   For example:

   .. jupyter-execute::

       import qse
       signal = qse.Signal([1, 1])
       signal = signal * 3 + 0.5
       print(signal)



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



   .. py:method:: time_per_value()

      Get the duration per value. Recall that the values are equally
      spaced over the total duration.

      :returns: *int* -- The duration of each value.



   .. py:method:: expand()

      Get an array of length 'duration' whose entries are the signal values.

      :returns: *np.ndarray* -- An array representing the signal.



   .. py:method:: to_pulser()

      Convert to a Pulser Waveform.

      :returns: *pulser.waveforms.Waveform* -- The waveform.



   .. py:method:: draw(time_units=None, signal_units=None, title=None)

      Draw the signal.

      :Parameters: * **time_units** (*str, optional*) -- The units of the duration.
                   * **signal_units** (*str, optional*) -- The units of the signal.
                   * **title** (*str, optional*) -- A title for the plot.



