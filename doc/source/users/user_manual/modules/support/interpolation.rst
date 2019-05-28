Linear interpolation
~~~~~~~~~~~~~~~~~~~~

Linear interpolation is implemented in GammaLib through the :doxy:`GNodeArray` 
class. This class contains a collection of nodes :math:`x_i` that may be
used to describe a functional relation :math:`y_i=f(x_i)`. The following 
code illustrates how the :doxy:`GNodeArray` class is used
(see ``examples/cpp/interpolate/interpolate.cpp`` for the source code):

**C++**

.. code-block:: cpp
   :linenos:

   double x_i[] = {1.0, 4.0, 6.0};
   double y_i[] = {8.0, 7.0, 2.0};
   GNodeArray nodes(3, x_i);
   for (double x = 0; x < 10.0; x += 0.5) {
       nodes.set_value(x);
       double y = y_i[nodes.inx_left()]  * nodes.wgt_left() + y_i[nodes.inx_right()] * nodes.wgt_right();
       std::cout << "x=" << x << " : y=" << y << std::endl;
   }

**Python**

.. code-block:: python
   :linenos:

   x_i   = [1.0, 4.0, 6.0]
   y_i   = [8.0, 7.0, 2.0]
   nodes = gammalib.GNodeArray(x_i)
   for x in range(20):
       nodes.set_value(x*0.5)
       y = y_i[nodes.inx_left()]  * nodes.wgt_left() + y_i[nodes.inx_right()] * nodes.wgt_right()
       print('x=%f : y=%f' % (x, y))


In line 1, the nodes :math:`x_i` at which the function values :math:`y_i`
are given are declared, the actual function values :math:`y_i` are
declared in line 2. In line 3, a node array is constructed from the
node values. Note that the actual function values are not part of the
node array, only the node values are in fact used by the :doxy:`GNodeArray`
class.

In lines 4-8, the function is interpolated at a number of values in the
interval :math:`[0,10[`. In line 5, the :math:`x` value is set at which
the interpolation should be done. The interpolation is then done in
line 6 using the formula

.. math::
   y = y_{i_{\rm left}} \times w_{i_{\rm left}} + y_{i_{\rm right}} \times w_{i_{\rm right}}

where :math:`i_{\rm left}` and :math:`i_{\rm right}` are the node indices
that encompass the :math:`x` value, and :math:`w_{i_{\rm left}}` and
:math:`w_{i_{\rm right}}` are the weights with which the function values 
:math:`y_{i_{\rm left}}` and :math:`y_{i_{\rm right}}` need to be multiplied
to obtain the interpolated value :math:`y`. Note that

.. math::
   w_{i_{\rm left}} + w_{i_{\rm right}} = 1

The method also works for extrapolation.
For :math:`x < x_0`, :math:`i_{\rm left}=0` and :math:`i_{\rm right}=1`,
while for :math:`x > x_{i_{\rm last}}`, :math:`i_{\rm left}=i_{\rm last}-1`
and :math:`i_{\rm right}=i_{\rm last}` (where :math:`i_{\rm last}` is the
index of the last node, which is :math:`2` in the example above).
The weights are set so that :math:`y` is extrapolated linearly.

It is obvious that :doxy:`GNodeArray` needs at least 2 node values to operate.
