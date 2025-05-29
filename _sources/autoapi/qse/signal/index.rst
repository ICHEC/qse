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

   Signal class represents a 1D signal with values and duration.

   The class supports arithmetic operations with other signals and
   scalars. Specifically:

   - Adding a scalar to a Signal returns a new Signal with modified
   values and the same duration. Say W is a signal, then W + 3.0 returns
   a signal with same duration, and values W.values + 3.0.

   - Adding two Signal instances concatenates their values and sums
   their durations. So if w1, w2 are instantiation of Signal, then
   w = w1 + w2 gives signal with concatenated values, i.e.,
   w.values = [w1.values, w2.values], and added duration.
   w.duration = w1.duration + w2.duration

   :ivar values: The values of the signal.
   :vartype values: ndarray
   :ivar duration: Duration of the signal.

   :vartype duration: int

   .. note::

      Currently, the object gets created for multi-dim arrays as well.
      However, it should be used for 1D only, we haven't made it useful
      or consistent for multi-dim usage.


   .. py:property:: duration
      :type: int


      Duration of the signal.

      :returns: The duration.
      :rtype: int


