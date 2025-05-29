qse.utils
=========

.. py:module:: qse.utils

.. autoapi-nested-parse::

   Utility functions for a variety of usage within and outside of QSE.



Functions
---------

.. autoapisummary::

   qse.utils.int2bin


Module Contents
---------------

.. py:function:: int2bin(x, width=32)

   converts an integer array to array of equivalent binary strings.
   Equivalent to:
       int2bin = np.vectorize(lambda x, width=16: np.binary_repr(x,width=width))
   However vectorize version is a bit slower compared to the one below.


