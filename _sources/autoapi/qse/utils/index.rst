qse.utils
=========

.. py:module:: qse.utils

.. autoapi-nested-parse::

   Utility functions for a variety of usage within and outside of QSE.



Functions
---------

.. autoapisummary::

   qse.utils.int2bin
   qse.utils.print_environment


Module Contents
---------------

.. py:function:: int2bin(x, width=32)

   Converts an integer array to array of equivalent binary strings.

   :Parameters: * **x** (*np.ndarray*) -- An array of integers.
                * **width** (*int*) -- The length of the binary strings.

   :returns: *np.ndarray* -- The array of binary strings.

   .. rubric:: Notes

   This function is equivalent to:

   >>> int2bin = np.vectorize(lambda x, width=16: np.binary_repr(x,width=width))

   However vectorize version is a bit slower compared to the one below.


.. py:function:: print_environment()

   Print the Python and qse version of the environment.


