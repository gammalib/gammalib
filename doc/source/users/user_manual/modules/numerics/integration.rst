Integration and derivatives
===========================

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
