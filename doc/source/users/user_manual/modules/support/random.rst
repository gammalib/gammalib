Random number generation
~~~~~~~~~~~~~~~~~~~~~~~~

Random number generation is widely used within GammaLib for drawing event 
distributions from functions.

**C++**

.. code-block:: cpp
   :linenos:

   GRan ran;                              // Construct random number generate instance
   double uniform     = ran.uniform();    // Uniform random number
   double exponential = ran.exp(3.7);     // Exponential random number
   double poisson     = ran.poisson(2.4); // Poisson random number
   double chi2        = ran.chisq2();     // Chi2 random number

**Python**

.. code-block:: python
   :linenos:

   ran         = gammalib.GRan()          # Construct random number generate instance
   uniform     = ran.uniform()            # Uniform random number
   exponential = ran.exp(3.7)             # Exponential random number
   poisson     = ran.poisson(2.4)         # Poisson random number
   chi2        = ran.chisq2()             # Chi2 random number

In line 1 a random number generator is allocated. If control over the
seed value of the random number generator is needed (for example to draw
different samples), you may specify the seed value upon construction:

**C++**

.. code-block:: cpp
   :linenos:

   unsigned long long int seed = 123456;
   GRan ran(seed);

**Python**

.. code-block:: python
   :linenos:

   seed = 123456
   ran  = gammalib.GRan(seed)

The :doxy:`GRan::uniform` method returns a random number between 0 and 1. The
:doxy:`GRan::exp` method returns a random number of the exponential law

.. math::
   p(x) = \lambda \exp( -\lambda x )

where :math:`\lambda` is the parameter of the distribution. In line 2
above, :math:`\lambda=3.7`. This method may be used to simulate the 
occurence time of an event, where :math:`\lambda` is the mean event rate.
Convsersely, :math:`1/\lambda` is the mean waiting time between events.

The :doxy:`GRan::poisson` method draws a random number from the Poisson 
distribution. You mya use this method to simulate the number of events
in case that a given mean number :math:`\lambda` of events is known.
In line 3 above, :math:`\lambda=2.4`.

The :doxy:`GRan::chisq2` method draws random numbers from the propability 
distribution

.. math::
   p(x) = \frac{1}{2\pi} x \exp \left( -\frac{1}{2} x^2 \right)

This method can be used to simulate the random radial offset of a measured
source position from the true source position, assuming an azimuthally
symmetric 2D Gaussian probability distribution.
