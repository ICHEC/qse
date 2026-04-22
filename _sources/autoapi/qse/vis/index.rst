qse.vis
=======

.. py:module:: qse.vis

.. autoapi-nested-parse::

   Visualize
   ---------

   This is the visualisation module of QSE package.
   It defines a number of helper functions to better
   visualise data.

   One of the most common function is the `draw` function
   to draw the qubits.



Functions
---------

.. autoapisummary::

   qse.vis.draw_qbits
   qse.vis.draw_amp_and_det
   qse.vis.draw_signal
   qse.vis.bar
   qse.vis.view_matrix


Package Contents
----------------

.. py:function:: draw_qbits(qbits, radius=None, show_labels=False, colouring=None, units=None, equal_aspect=True)

   Visualize the positions of a set of qubits.

   :Parameters: * **qbits** (*qse.Qbits*) -- The Qbits object.
                * **radius** (*float | str*) -- A cutoff radius for visualizing bonds.
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


.. py:function:: draw_amp_and_det(amplitude, detuning, time_units=None, signal_units=None, title=None)

   Draw a signal amplitude together with a signal detuning.

   :Parameters: * **amplitude** (*qse.Signal | qse.Signals*) -- The amplitude to be drawn.
                * **detuning** (*qse.Signal | qse.Signals*) -- The detuning to be drawn.
                * **time_units** (*str, optional*) -- The units of the duration.
                * **signal_units** (*str, optional*) -- The units of the signal.
                * **title** (*str, optional*) -- A title for the plot.


.. py:function:: draw_signal(signal, time_units=None, signal_units=None, title=None)

   Draw the signal.

   :Parameters: * **signal** (*qse.Signal | qse.Signals*) -- The signal to be drawn.
                * **time_units** (*str, optional*) -- The units of the duration.
                * **signal_units** (*str, optional*) -- The units of the signal.
                * **title** (*str, optional*) -- A title for the plot.


.. py:function:: bar(dict, cutoff=0, ylabel='Count')

   Plot a bar chart from a dictionary, filtering values below a cutoff.

   :Parameters: * **dict** (*dict*) -- A dictionary where keys are categories and values are their corresponding
                  counts or values.
                * **cutoff** (*int or float, optional*) -- Minimum value threshold for inclusion in the plot. Keys with values less
                  than or equal to `cutoff` are excluded.
                  Default is 0.
                * **ylabel** (*str, optional*) -- Label for the y-axis of the plot. Default is "Count".


.. py:function:: view_matrix(matrix, labels_x=None, labels_y=None, vcenter=None)

   Visualise a matrix.

   :Parameters: * **matrix** (*np.ndarray*) -- The matrix to be visualised.
                * **labels_x** (*list, optional*) -- Labels to be displayed on the x axis.
                * **labels_y** (*list, optional*) -- Labels to be displayed on the y axis.
                * **vcenter** (*float, optional*) -- The center of the colorbar.


