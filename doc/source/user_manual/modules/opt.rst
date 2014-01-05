.. _sec_opt:

Optimizers
----------

Overview
~~~~~~~~

The optimizer module provides classes for function optimization.
The abstract ``GOptimizerFunction`` base class defines the interface
for the function that should be optimized. The ``GOptimizerPar``
class defines a parameter of the function, and the ``GOptimizerPars``
container class collects all function parameters. The optimizer
is represented by the abstract ``GOptimizer`` base class. So far,
the only optimizer algorithm implemented in GammaLib is the
Levenberg-Marquardt algorithm, implemented by the class
``GOptimizerLM``.


:ref:`fig_uml_opt` presents an overview over the C++ classes of
the optimizer module and their relations.

.. _fig_uml_opt:

.. figure:: uml_opt.png
   :width: 100%

   Optimizer module

The model parameter class ``GModelPar``, as a specific implementation
of the optimizer parameter class, derives from ``GOptimizerPar``.
One implementation of an optimizer function is the
``GObservations::likelihood`` class that is used for maximum
likelihood fitting within GammaLib.


Optimizing a function
~~~~~~~~~~~~~~~~~~~~~

This example illustrates how the minimum of a quadratic function of 
the form :math:`f(x)=2ax^2+bx+c` can be determined using the optimizer.
Please see ``examples/cpp/optimize`` for the source code.

First, the function is implemented as a class derived from
``GOptimizerFunction``::

    1  class function : public GOptimizerFunction {
    2  public:
    3      function(void) : m_value(0), m_gradient(1), m_covar(1,1) {}
    4      void           eval(const GOptimizerPars& pars);
    5      double         value(void) { return m_value; }
    6      GVector*       gradient(void) { return &m_gradient; }
    7      GMatrixSparse* covar(void) { return &m_covar; }
    8  protected:
    9      double        m_value;    //!< Function value
   10      GVector       m_gradient; //!< Function gradient vector
   11      GMatrixSparse m_covar;    //!< Covariance matrix
   12  };
   13  void function::eval(const GOptimizerPars& pars)
   14  {
   15      const double a =  2.0;
   16      const double b = -4.0;
   17      const double c =  2.0;
   18      double x       = pars[0]->value();
   19      m_value        = a*x*x + b*x + c;
   20      m_gradient[0]  = 2.0*a*x + b;
   21      m_covar(0,0)   = m_gradient[0] * m_gradient[0];
   22  };

Lines 1-12 define the ``function`` class that requires implementing the
``eval``, ``value``, ``gradient`` and ``covar`` methods. The ``eval``
method, implement in lines 13-22, performs the computation of the
function value, the gradient and the products of the gradients.

The optimization is then done using the following code::

    1  GOptimizerLM opt;
    2  function fct;
    3  GOptimizerPars pars(1);
    4  pars[0]->value(1.5);
    5  opt.optimize(fct, pars);
    6  std::cout << "Function value .....: " << fct.value() << std::endl;
    7  std::cout << "Parameter value ....: " << pars[0]->value() << std::endl;

Line 1 allocates an Levenberg-Marquardt optimizer, line 2 creates an 
instance of the function to optimize. In line 3, a parameter container
with a single parameter is allocated, and the value of the single parameter
is set to 1.5 in line 4. In line 5, the optimizer is called, and the 
resulting function value and best fitted parameter is logged in to console 
in lines 6-7.


