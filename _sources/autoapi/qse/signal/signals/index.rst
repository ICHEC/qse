qse.signal.signals
==================

.. py:module:: qse.signal.signals

.. autoapi-nested-parse::

   Definition of the Signal class.



Classes
-------

.. autoapisummary::

   qse.signal.signals.Signals


Module Contents
---------------

.. py:class:: Signals(signals=None)

   The Signal class represents a collection of Signals.

   :Parameters: **signals** (*list[qse.Signal], optional*) -- The signals.

   .. rubric:: Examples

   .. jupyter-execute::

       import qse
       import numpy as np

       x = qse.Signal([1], 10)
       y = qse.Signal(np.linspace(0, 1, 5), 5)
       ss = qse.Signals([x, y])
       print(ss)
       z = qse.Signals()
       z += x
       z += y
       print(z)

   As shown above, one can also create Signals by adding two signals.


   .. py:method:: to_pulser()

      Convert to a Pulser Waveform.

      :returns: *pulser.waveforms.Waveform* -- The waveform.



   .. py:method:: expand()

      Get an array of length 'duration' whose entries are the signal values.

      :returns: *np.ndarray* -- An array representing the signal.



   .. py:method:: draw(time_units=None, signal_units=None, title=None)

      Draw the signal.

      :Parameters: * **time_units** (*str, optional*) -- The units of the duration.
                   * **signal_units** (*str, optional*) -- The units of the signal.
                   * **title** (*str, optional*) -- A title for the plot.



