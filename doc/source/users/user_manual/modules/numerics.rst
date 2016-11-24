.. _sec_numerics:

Numerical methods
-----------------

Overview
~~~~~~~~

This module contains classes and functions that are needed for numerical 
computations within GammaLib. The module provides support for differentiation
and integration of one-dimensional functions. The interface for the 
one-dimensional functions is defined by the abstract :doxy:`GFunction` base
class, integration is performed using the :doxy:`GIntegral` class, 
differentiation is done using the :doxy:`GDerivative` class. In addition,
numerical constants and function that are extensively used throughout GammaLib
are defined in the :doxy:`GMath.hpp` header file. All constants and functions
are declared in the :doxy:`gammalib` namespace.


Constants
~~~~~~~~~

The following constants are available:

========================== =====
Constant                   Value
========================== =====
``gammalib::pi``           :math:`\pi`
``gammalib::twopi``        :math:`2\pi`
``gammalib::fourpi``       :math:`4\pi`
``gammalib::pihalf``       :math:`\pi/2`
``gammalib::inv_pihalf``   :math:`(\pi/2)^{-1}`
``gammalib::inv_sqrt4pi``  :math:`(\sqrt{4\pi})^{-1}`
``gammalib::pi2``          :math:`\pi^2`
``gammalib::deg2rad``      :math:`\pi/180` (multiplication converts degrees to radians)
``gammalib::rad2deg``      :math:`180/\pi` (multiplication converts radians to degrees)
``gammalib::ln2``          :math:`\log 2` (natural logarithm of 2)
``gammalib::ln10``         :math:`\log 10` (natural logarithm of 10)
``gammalib::inv_ln2``      :math:`(\log 2)^{-1}`
``gammalib::inv_ln10``     :math:`(\log 10)^{-1}`
``gammalib::onethird``     :math:`1/3`
``gammalib::twothird``     :math:`2/3`
``gammalib::fourthird``    :math:`4/3`
``gammalib::sqrt_onehalf`` :math:`\sqrt{1/2}`
``gammalib::sqrt_pihalf``  :math:`\sqrt{\pi/2}`;
``gammalib::sqrt_two``     :math:`\sqrt{2}`;
========================== =====


Functions
~~~~~~~~~

The following functions are available:

===================== ===========
Function              Description
===================== ===========
``gammalib::acos``    Arc cosine function that avoids NaN due to rounding errors
``gammalib::cosd``    Cosine function for argument given in degrees
``gammalib::sind``    Sine function for argument given in degrees
``gammalib::tand``    Tangens function for argument given in degrees
``gammalib::asind``   Arc sine function for argument given in degrees
``gammalib::acosd``   Arc cosine function for argument given in degrees
``gammalib::atand``   Arc tangens function for argument given in degrees
``gammalib::atan2d``  Arc tangens function for argument x/y given in degrees
``gammalib::sincosd`` Returns sine and cosine for argument given in degrees
``gammalib::gammln``  Natural logarithm of gamma function
``gammalib::erfcc``   Complementary error function
``gammalib::erfinv``  Inverse error function
``gammalib::modulo``  Remainder of division x/y
===================== ===========


Integration and derivatives
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following code illustrates how integrations and derivatives are
computed within GammaLib (see ``examples/cpp/numerics/numerics.cpp`` for the
source code):

.. code-block:: cpp
   :linenos:

    class function : public GFunction {
    public:
        function(const double& sigma) : m_a(1.0/(sigma*std::sqrt(gammalib::twopi))), m_sigma(sigma) {}
        double eval(const double& x) { return m_a*std::exp(-0.5*x*x/(m_sigma*m_sigma)); }
    protected:
        double m_a;
        double m_sigma;
    };
    int main(void) {
        function fct(3.0);
        GIntegral integral(&fct);
        integral.eps(1.0e-8);
        double result = integral.romb(-15.0, +15.0);
        std::cout << "Integral:       " << result << std::endl;
        GDerivative derivative(&fct);
        std::cout << "Derivative(0):  " << derivative.value(0.0) << std::endl;
    return 0;
    }

The function that should be integrated or differentiated is defined in
lines 1-8 as a class that derives from the abstract :doxy:`GFunction` base
class. The only method that needs to be implement in the derived class,
here named ``function`` is the :doxy:`GFunction::eval` method that takes a const reference
to a double precision value as argument and that returns a double precision
value, which is the function value evaluated at the argument. Parameters
may be passed to the function upon construction, as illustrated by the
``m_a`` and ``m_sigma`` members that are initialised by the constructor.

The function is allocated in line 10 with a sigma parameter of 3. Line 11
the prepares for the integration by allocating an integration object. The
:doxy:`GIntegral` constructor takes a reference to the function as argument.
In line 12, the relative precision of the integration object is set to
:math:`10^{-8}` (by default the precision is set to :math:`10^{-6}`).
In line 13, the integration is done over the parameter interval
:math:`[-15,15]`. As this covers basically the entire area of the
Gaussian function, the result will be very close to 1 (the result is
printed in line 14). Note that the Romberg method is used for integration
by invoking the ``romb`` method. This is the only method that is so far
available in GammaLib.

Differentiating a function is similar. For this purpose, a :doxy:`GDerivative`
object is created in line 15 with takes a reference to the function as
argument. Using the :doxy:`GDerivative::value` method, the derivative is computed in line
16 for a function argument of 0. As the Gaussian has a maximum there, the
result will be 0.
