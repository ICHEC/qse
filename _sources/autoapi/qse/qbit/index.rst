qse.qbit
========

.. py:module:: qse.qbit

.. autoapi-nested-parse::

   Functions adapted from ASE's Qbit/Qbits styled objects.



Classes
-------

.. autoapisummary::

   qse.qbit.Qbit


Functions
---------

.. autoapisummary::

   qse.qbit.qbitproperty
   qse.qbit.xyzproperty


Module Contents
---------------

.. py:function:: qbitproperty(name, doc)

   Helper function to easily create Qbit attribute property.


.. py:function:: xyzproperty(index)

   Helper function to easily create Qbit XYZ-property.


.. py:class:: Qbit(label='X', state=(1, 0), position=(0, 0, 0), qbits=None, index=None)

   Class for representing a single qbit.

   :Parameters: * **label** (*str or int*) -- Can be a str or an int label.
                * **state** (*list or tuple or np.ndarray*) -- Quantum state of the qubit
                * **position** (*list or np.ndarray*) -- Sequence of 3 floats qubit position.
                * **qbits** (*qse.Qbits*) -- The Qbits object that the Qbit is attached
                  to. Defaults to None.
                * **index** (*int*) -- The associated index of the Qbit in the
                  Qbits object. Defaults to None.

   .. rubric:: Notes

   You can create a qbit object with

   >>> q = qse.Qbit()


   .. py:method:: cut_reference_to_qbits()

      Cut reference to qbits object.



   .. py:method:: get_raw(name)

      Get name attribute, return None if not explicitly set.



   .. py:method:: get(name)

      Get name attribute, return default if not explicitly set.



   .. py:method:: set(name, value)

      Set name attribute to value.



   .. py:method:: delete(name)

      Delete name attribute.



