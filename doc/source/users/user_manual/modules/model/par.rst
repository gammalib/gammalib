Model parameters
================

Model parameters are implemented by :doxy:`GModelPar` which derives from
the abstract :doxy:`GOptimizerPar` base class. Each model parameter is
factorised into a value and a scale part according to

.. math::
   p = v \times s

where
:math:`p` is the model parameter,
:math:`v` is the model value factor, and
:math:`s` is the model scale factor.
All three components can be accessed or set via the
:doxy:`GModelPar::value()`,
:doxy:`GModelPar::factor_value()`, and
:doxy:`GModelPar::scale()` methods.

In addition, the class handles parameter errors and gradients as well
as parameter minima or maxima.
Parameter errors can be accessed or set via the
:doxy:`GModelPar::error()` and :doxy:`GModelPar::factor_error()`,
where the former concerns the error on :math:`p` while the latter concerns
the error on :math:`v`.
Equivalent methods exist for
the gradient
(:doxy:`GModelPar::gradient()` and :doxy:`GModelPar::factor_gradient()`),
the minimum
(:doxy:`GModelPar::min()` and :doxy:`GModelPar::factor_min()`), and
the maximum
(:doxy:`GModelPar::max()` and :doxy:`GModelPar::factor_max()`).

Finally, :doxy:`GModelPar::fix()` and :doxy:`GModelPar::free()` methods
allow to fix or to free a model parameter, and the :doxy:`GModelPar::is_free()`
and :doxy:`GModelPar::is_fix()` methods allow to check whether a parameter is
free or fixed.


