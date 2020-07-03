.. _um_cta_background:

Modelling CTA background
========================

The standard method for modelling the CTA instrumental background makes 
use of background rate information that is included in the instrument 
response functions.
For unbinned or binned analysis the background model is implemented by
the :doxy:`GCTAModelIrfBackground` class that derives from the
abstract :doxy:`GModelData` class.
The XML description of the background model has the following format:

.. code-block:: xml

   <source name="Background" type="CTAIrfBackground" instrument="CTA">
     <spectrum type="PowerLaw">
       <parameter name="Prefactor" scale="1.0" value="1.0" min="1e-3" max="1e3"    free="1"/>
       <parameter name="Index"     scale="1.0" value="0.0" min="-5.0" max="+5.0"   free="1"/>
       <parameter name="Scale"     scale="1e6" value="1.0" min="0.01" max="1000.0" free="0"/>
     </spectrum>
   </source>

The type of the model is ``CTAIrfBackground``.
The tag ``<spectrum>`` allows to specify a spectral model with which the
background model in the IRF will be multiplied.
Any spectral model available in GammaLib can be used as the spectral
component of the model.
This allows to adjust the energy distribution of the background template
to the data.
In the above example the background rates are multiplied with a power law.

.. note::

  In case that no background information is available in the instrument
  response function a first order background model can be obtained by
  using the effective area as a spatial template for the background.
  The XML model description has the following format:

  .. code-block:: xml
  
     <source name="Background" type="CTAAeffBackground" instrument="CTA">
       <spectrum type="PowerLaw">
         <parameter name="Prefactor" scale="1e-14" value="1.0"  min="1e-3" max="1e3"    free="1"/>
         <parameter name="Index"     scale="1.0"   value="-2.4" min="-5.0" max="+5.0"   free="1"/>
         <parameter name="Scale"     scale="1e6"   value="1.0"  min="0.01" max="1000.0" free="0"/>
       </spectrum>
     </source>

For stacked analysis there is an equivalent model that makes use of the 
background model provided in the ``BkgCube`` parameter of the observation
definition XML file:

.. code-block:: xml

   <source name="Background" type="CTACubeBackground" instrument="CTA">
     <spectrum type="PowerLaw">
       <parameter name="Prefactor" scale="1.0" value="1.0" min="1e-3" max="1e3"    free="1"/>
       <parameter name="Index"     scale="1.0" value="0.0" min="-5.0" max="+5.0"   free="1"/>
       <parameter name="Scale"     scale="1e6" value="1.0" min="0.01" max="1000.0" free="0"/>
     </spectrum>
   </source>

The structure of the XML model is identical to the one shown before, but 
the type of the model is now ``CTACubeBackground``.

Another type of background model exists that assumes that the 
spatial and spectral component of the background can be factorised and 
that the spatial component is radially symmetric.
This background model is rather simplistic, and is mainly still there for 
historic reasons.
Here an example of the XML definition.
The type of the model is ``RadialAcceptance``.

.. code-block:: xml

   <source name="Background" type="RadialAcceptance" instrument="CTA">
     <spectrum type="PowerLaw">
       <parameter name="Prefactor" scale="1e-6" value="61.8" min="0.0"  max="1000.0" free="1"/>
       <parameter name="Index"     scale="-1"   value="1.85" min="0.0"  max="+5.0"   free="1"/>
       <parameter name="Scale"     scale="1e6"  value="1.0"  min="0.01" max="1000.0" free="0"/>
     </spectrum>
     <radialModel type="Gaussian">
       <parameter name="Sigma" scale="1.0" value="3.0" min="0.01" max="10.0" free="1"/>
     </radialModel>
   </source>

The spatial component of the model is a Gaussian in offset angle squared.
Alternativly, a profile can be specified:

.. code-block:: xml

   <source name="Background" type="RadialAcceptance" instrument="CTA">
     <spectrum type="PowerLaw">
       <parameter name="Prefactor" scale="1e-6" value="61.8" min="0.0"  max="1000.0" free="1"/>
       <parameter name="Index"     scale="-1"   value="1.85" min="0.0"  max="+5.0"   free="1"/>
       <parameter name="Scale"     scale="1e6"  value="1.0"  min="0.01" max="1000.0" free="0"/>
     </spectrum>
     <radialModel type="Profile">
       <parameter name="Width" scale="1.0" value="1.5" min="0.1" max="1000.0" free="1"/>
       <parameter name="Core"  scale="1.0" value="3.0" min="0.1" max="1000.0" free="1"/>
       <parameter name="Tail"  scale="1.0" value="5.0" min="0.1" max="1000.0" free="1"/>
     </radialModel>
   </source>

Or a polynom:

.. code-block:: xml

   <source name="Background" type="RadialAcceptance" instrument="CTA">
     <spectrum type="PowerLaw">
       <parameter name="Prefactor" scale="1e-6" value="61.8" min="0.0"  max="1000.0" free="1"/>
       <parameter name="Index"     scale="-1"   value="1.85" min="0.0"  max="+5.0"   free="1"/>
       <parameter name="Scale"     scale="1e6"  value="1.0"  min="0.01" max="1000.0" free="0"/>
     </spectrum>
     <radialModel type="Polynom">
       <parameter name="Coeff0" scale="1.0" value="+1.00000"   min="-10.0" max="10.0" free="0"/>
       <parameter name="Coeff1" scale="1.0" value="-0.1239176" min="-10.0" max="10.0" free="1"/>
       <parameter name="Coeff2" scale="1.0" value="+0.9751791" min="-10.0" max="10.0" free="1"/>
       <parameter name="Coeff3" scale="1.0" value="-3.0584577" min="-10.0" max="10.0" free="1"/>
       <parameter name="Coeff4" scale="1.0" value="+2.9089535" min="-10.0" max="10.0" free="1"/>
       <parameter name="Coeff5" scale="1.0" value="-1.3535372" min="-10.0" max="10.0" free="1"/>
       <parameter name="Coeff6" scale="1.0" value="+0.3413752" min="-10.0" max="10.0" free="1"/>
       <parameter name="Coeff7" scale="1.0" value="-0.0449642" min="-10.0" max="10.0" free="1"/>
       <parameter name="Coeff8" scale="1.0" value="+0.0024321" min="-10.0" max="10.0" free="1"/>
     </radialModel>
   </source>
