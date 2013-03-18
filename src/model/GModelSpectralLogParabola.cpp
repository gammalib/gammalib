/***************************************************************************
 *    GModelSpectralLogParabola.cpp - Log parabola spectral model class    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2013 by Michael Mayer                               *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *  This program is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU General Public License as published by   *
 *  the Free Software Foundation, either version 3 of the License, or      *
 *  (at your option) any later version.                                    *
 *                                                                         *
 *  This program is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *  GNU General Public License for more details.                           *
 *                                                                         *
 *  You should have received a copy of the GNU General Public License      *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GModelSpectralLogParabola.cpp
 * @brief Log parabola spectral model class definition
 * @author Michael Mayer 
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include "GException.hpp"
#include "GTools.hpp"
#include "GIntegral.hpp"
#include "GModelSpectralLogParabola.hpp"
#include "GModelSpectralRegistry.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GModelSpectralLogParabola g_spectral_logparabola_seed;
const GModelSpectralRegistry    g_spectral_logparabola_registry(&g_spectral_logparabola_seed);

/* __ Method name definitions ____________________________________________ */
#define G_FLUX          "GModelSpectralLogParabola::flux(GEnergy&, GEnergy&)"
#define G_EFLUX        "GModelSpectralLogParabola::eflux(GEnergy&, GEnergy&)"
#define G_MC       "GModelSpectralLogParabola::mc(GEnergy&, GEnergy&, GRan&)"
#define G_READ                "GModelSpectralLogParabola::read(GXmlElement&)"
#define G_WRITE              "GModelSpectralLogParabola::write(GXmlElement&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GModelSpectralLogParabola::GModelSpectralLogParabola(void) : GModelSpectral()
{
    // Initialise private members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor
 *
 * @param[in] norm Power law normalization.
 * @param[in] index Power law index.
 * @param[in] curvature Curvature.
 *
 * Construct a LogParabola model from a normalization value, a spectral
 * index and a curvature
 ***************************************************************************/
GModelSpectralLogParabola::GModelSpectralLogParabola(const double& norm,
                                                     const double& index,
                                                     const double& curvature) :
                           GModelSpectral()
{
    // Initialise members
    init_members();

    // Set parameters
    m_norm.Value(norm);
    m_index.Value(index);
    m_curvature.Value(curvature);

    // Autoscale parameters
    autoscale();

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 *
 * Creates instance of a log parabola spectral model by extracting information
 * from an XML element. See GModelSpectralLogParabola::read() for more information
 * about the expected structure of the XML element.
 ***************************************************************************/
GModelSpectralLogParabola::GModelSpectralLogParabola(const GXmlElement& xml) : 
                           GModelSpectral()
{
    // Initialise members
    init_members();

    // Read information from XML element
    read(xml);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] model LogParabola model.
 ***************************************************************************/
GModelSpectralLogParabola::GModelSpectralLogParabola(const GModelSpectralLogParabola& model) :
                           GModelSpectral(model)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(model);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GModelSpectralLogParabola::~GModelSpectralLogParabola(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Operators                                =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] model LogParabola model.
 * @return LogParabola model.
 ***************************************************************************/
GModelSpectralLogParabola& GModelSpectralLogParabola::operator=(const GModelSpectralLogParabola& model)
{
    // Execute only if object is not identical
    if (this != &model) {

        // Copy base class members
        this->GModelSpectral::operator=(model);

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(model);

    } // endif: object was not identical

    // Return
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                              Public methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear instance
 ***************************************************************************/
void GModelSpectralLogParabola::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GModelSpectral::free_members();

    // Initialise members
    this->GModelSpectral::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
 *
 * @return Pointer to deep copy of log parabola spectrum.
 ***************************************************************************/
GModelSpectralLogParabola* GModelSpectralLogParabola::clone(void) const
{
    return new GModelSpectralLogParabola(*this);
}


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] srcEng True energy of photon.
 * @return Function value.
 *
 * The logparabola function is defined as
 * \f[I(E)=norm (E/pivot)^{index-curvature*log(E/pivot)}\f]
 * where
 * \f$norm\f$ is the normalization or prefactor,
 * \f$pivot\f$ is the pivot energy, and
 * \f$index\f$ is the spectral index.
 * \f$curvature\f$ is the curvature
 *
 * @todo For the moment the pivot energy is fixed to units of MeV. This may
 * not be ideal and should eventually be improved in the future.
 * Furthermore, the method expects that energy!=0. Otherwise Inf or NaN
 * may result.
 ***************************************************************************/
double GModelSpectralLogParabola::eval(const GEnergy& srcEng) const
{
    // Compute function value
    double energy   = srcEng.MeV() / pivot();
    double exponent = index() + curvature()*log(energy);
    double power    = std::pow(energy,exponent);
    double value    = norm() * power;

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (isnotanumber(value) || isinfinite(value)) {
        std::cout << "*** ERROR: GModelSpectralLogParabola::eval";
        std::cout << "(srcEng=" << srcEng << "):";
        std::cout << " NaN/Inf encountered";
        std::cout << " (value=" << value;
        std::cout << ", energy=" << energy;
        std::cout << ", index=" << index();
        std::cout << ", curvature=" << curvature();
        std::cout << ", pivot=" << pivot();
        std::cout << ", power=" << power;
        std::cout << ")" << std::endl;
    }
    #endif

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Evaluate function and gradients
 *
 * @param[in] srcEng True energy of photon.
 * @return Function value.
 *
 * The logparabol function is defined as
 * \f[I(E)=norm (E/pivot)^{index-curvature\ln(E/pivot)}\f]
 * where
 * \f$norm=n_s n_v\f$ is the normalization or prefactor,
 * \f$pivot=p_s p_v\f$ is the pivot energy, and
 * \f$index=i_s i_v\f$ is the spectral index.
 * \f$curvature=i_s i_v\f$ is the curvature
 * Note that each parameter is factorised into a scaling factor and a value
 * and that the method is expected to return the gradient with respect to
 * the parameter value (i.e. n_v, p_v, i_v, and c_v in this case).
 *
 * The partial derivatives of the parameter values are given by
 *@todo implement formulae here
 *
 * @todo For the moment the pivot energy is fixed to units of MeV. This may
 * not be ideal and should eventually be improved in the futur.
 * Furthermore, the method expects that energy!=0. Otherwise Inf or NaN
 * may result.
 ***************************************************************************/
double GModelSpectralLogParabola::eval_gradients(const GEnergy& srcEng) const
{
    // Compute function value
    double energy   = srcEng.MeV() / pivot();
    double exponent = index() + curvature()*std::log(energy);
    double power    = std::pow(energy,exponent);
    double value    = norm() * power;

    // Compute partial derivatives of the parameter values
    double log_energy  = std::log(energy);
    
    double g_norm      = (m_norm.isfree())  ? m_norm.scale() * power : 0.0;
    double g_index     = (m_index.isfree())
                         ? value * m_index.scale() * log_energy : 0.0;
    double g_curvature = (m_curvature.isfree())
                         ? value * m_curvature.scale() * log_energy * 
                           log_energy
                         : 0.0;
    double g_pivot     = (m_pivot.isfree())
                         ? -value/m_pivot.factor_value() * (exponent + curvature() *
                           log_energy)
                         : 0.0;

    // Set gradients (circumvent const correctness)
    const_cast<GModelSpectralLogParabola*>(this)->m_norm.factor_gradient(g_norm);
    const_cast<GModelSpectralLogParabola*>(this)->m_index.factor_gradient(g_index);
    const_cast<GModelSpectralLogParabola*>(this)->m_curvature.factor_gradient(g_curvature);
    const_cast<GModelSpectralLogParabola*>(this)->m_pivot.factor_gradient(g_pivot);

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (isnotanumber(value) || isinfinite(value)) {
        std::cout << "*** ERROR: GModelSpectralLogParabola::eval_gradients";
        std::cout << "(srcEng=" << srcEng << "):";
        std::cout << " NaN/Inf encountered";
        std::cout << " (value=" << value;
        std::cout << ", energy=" << energy;
        std::cout << ", index=" << index();
        std::cout << ", curvature=" << curvature();
        std::cout << ", pivot=" << pivot();
        std::cout << ", power=" << power;
        std::cout << ", g_norm=" << g_norm;
        std::cout << ", g_index=" << g_index;
        std::cout << ", g_curvature=" << g_curvature;
        std::cout << ", g_pivot=" << g_pivot;
        std::cout << ")" << std::endl;
    }
    #endif

    // Return
    return value;
}


/***********************************************************************//**
 * @brief Returns model photon flux between [emin, emax] (units: ph/cm2/s)
 *
 * @param[in] emin Minimum photon energy.
 * @param[in] emax Maximum photon energy.
 * @return Photon flux (ph/cm2/s).
 *
 * Computes
 * \f[\int_{E_{\rm min}}^{E_{\rm max}} I(E) dE\f]
 * where
 * \f$E_{\rm min}\f$ and \f$E_{\rm max}\f$ are the minimum and maximum
 * energy, respectively, and
 * \f$I(E)\f$ is the spectral model (units: ph/cm2/s/MeV).
 * The integration is done numerically.
 ***************************************************************************/
double GModelSpectralLogParabola::flux(const GEnergy& emin,
                                       const GEnergy& emax) const
{
	// Initialise function to integrate
    flux_kern flux(norm(),index(),curvature(),pivot());

    // Initialise integral class with function
    GIntegral integral(&flux);

    // Set integration precision
    integral.eps(1.0e-8);

    // Calculate integral between emin and emax
    double photonflux = integral.romb(emin.MeV(), emax.MeV());

    //Return value
    return photonflux;
}


/***********************************************************************//**
 * @brief Returns model energy flux between [emin, emax] (units: erg/cm2/s)
 *
 * @param[in] emin Minimum photon energy.
 * @param[in] emax Maximum photon energy.
 * @return Energy flux (erg/cm2/s).
 *
 * Computes
 * \f[\int_{E_{\rm min}}^{E_{\rm max}} I(E) E dE\f]
 * where
 * \f$E_{\rm min}\f$ and \f$E_{\rm max}\f$ are the minimum and maximum
 * energy, respectively, and
 * \f$I(E)\f$ is the spectral model (units: ph/cm2/s/MeV).
 * The integration is done numerically.
 ***************************************************************************/
double GModelSpectralLogParabola::eflux(const GEnergy& emin,
                                        const GEnergy& emax) const
{
	// Initialise function to integrate
    eflux_kern eflux(norm(),index(),curvature(),pivot());

    // Initialise integral class with function
    GIntegral integral(&eflux);

    // Set integration precision
    integral.eps(1.0e-8);

    // Calculate integral between emin and emax
    double energyflux = integral.romb(emin.MeV(), emax.MeV());

    // Return value
    return energyflux;
}


/***********************************************************************//**
 * @brief Returns Monte Carlo energy between [emin, emax]
 *
 * @param[in] emin Minimum photon energy.
 * @param[in] emax Maximum photon energy.
 * @param[in] ran Random number generator.
 * @return Energy.
 *
 * @exception GException::feature_not_implemented
 *            Feature not yet implemented.
 *
 * Returns Monte Carlo energy by randomly drawing from a power law.
 *
 ***************************************************************************/
GEnergy GModelSpectralLogParabola::mc(const GEnergy& emin,
                                      const GEnergy& emax,
                                      GRan& ran) const
{
		if (emin >= emax) {
	        throw GException::erange_invalid(G_MC, emin.MeV(), emax.MeV(),
	              "Minimum energy < maximum energy required.");
	    }

	    // Allocate energy
	    GEnergy energy;

	    // Update cache
	    update_mc_cache(emin, emax);

	    // Initialise energy
	    double eng;

	    // Initialse acceptance fraction
	    double acceptance_fraction;

	    // Use rejection method to draw a random energy. We first draw
	    // analytically from a power law, and then compare the power law
	    // at the drawn energy to the curved function. This
	    // gives an acceptance fraction, and we accept the energy only if
	    // a uniform random number is <= the acceptance fraction.
	    do {

	        // Get uniform random number
	        double u = ran.uniform();

	        // Case A: Corresponding mc Plaw-Index is not -1
	        if (m_mc_exponent != 0.0) {
	            if (u > 0.0) {
	                eng = std::exp(std::log(u * m_mc_pow_ewidth + m_mc_pow_emin) /
	                               m_mc_exponent);
	            }
	            else {
	                eng = 0.0;
	            }
	        }

	        // Case B: Corresponding mc Plaw-Index is  -1
	        else {
	            eng = std::exp(u * m_mc_pow_ewidth + m_mc_pow_emin);
	        }

	        // Compute powerlaw at given energy
	        double e_norm = eng / pivot();
	        double plaw   = m_mc_norm*std::pow(eng/pivot(), m_mc_exponent-1.0);

	        // Compute logparabola at given energy
	        double logparabola = norm()*std::pow(e_norm,index()+curvature()*std::log(e_norm));

	        // Compute acceptance fraction
	        acceptance_fraction = logparabola / plaw;

	    } while (ran.uniform() > acceptance_fraction);

	    // Set energy
	    energy.MeV(eng);

	    // Return energy
	    return energy;
}


/***********************************************************************//**
 * @brief Read model from XML element
 *
 * @param[in] xml XML element containing logparabola model information.
 *
 * @exception GException::model_invalid_parnum
 *            Invalid number of model parameters found in XML element.
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter names found in XML element.
 *
 * Read the spectral log parabola information from an XML element. The format
 * of the XML elements is
 *
 *     <spectrum type="LogParabola">
 *       <parameter name="Prefactor" scale=".." value=".." min=".." max=".." free=".."/>
 *       <parameter name="Index"     scale=".." value=".." min=".." max=".." free=".."/>
 *       <parameter name="Curvature" scale=".." value=".." min=".." max=".." free=".."/>
 *       <parameter name="Scale"     scale=".." value=".." min=".." max=".." free=".."/>
 *     </spectrum>
 *
 * or for compliance with Fermi-LAT
 *
 *     <spectrum type="LogParabola">
 *       <parameter name="norm"  scale=".." value=".." min=".." max=".." free=".."/>
 *       <parameter name="alpha" scale=".." value=".." min=".." max=".." free=".."/>
 *       <parameter name="beta"  scale=".." value=".." min=".." max=".." free=".."/>
 *       <parameter name="Eb"    scale=".." value=".." min=".." max=".." free=".."/>
 *     </spectrum>
 *
 * @todo Add parameter validity check
 ***************************************************************************/
void GModelSpectralLogParabola::read(const GXmlElement& xml)
{
    // Verify that XML element has exactly 3 parameters
    if (xml.elements() != 4 || xml.elements("parameter") != 4) {
        throw GException::model_invalid_parnum(G_READ, xml,
              "LogParabola model requires exactly 4 parameters.");
    }

    // Extract model parameters
    int npar[] = {0, 0, 0,0};
    for (int i = 0; i < 4; ++i) {

        // Get parameter element
        const GXmlElement* par = xml.element("parameter", i);

        // Handle prefactor
        if (par->attribute("name") == "Prefactor" ||
        	par->attribute("name") == "norm") {
            m_norm.read(*par);
            npar[0]++;
        }

        // Handle index
        else if (par->attribute("name") == "Index"){
            m_index.read(*par);
            npar[1]++;
        }

        // Change sign if index is defined Fermi-like
        else if(par->attribute("name") == "alpha") {
        	m_index.read(*par);
        	m_index.Value(-m_index.Value());
        	npar[1]++;
        }

        // Handle curvature
        else if (par->attribute("name") == "Curvature") {
            m_curvature.read(*par);
            npar[2]++;
        }

        // Change sign if curvature is defined Fermi-like
        else if(par->attribute("name") == "beta") {
        	m_curvature.read(*par);
        	m_curvature.Value(-m_curvature.Value());
        	npar[2]++;
        }

        // Handle pivot energy
        else if (par->attribute("name") == "Scale" ||
			     par->attribute("name") == "Eb") {
            m_pivot.read(*par);
            npar[3]++;
        }

    } // endfor: looped over all parameters

    // Verify that all parameters were found
    if (npar[0] != 1 || npar[1] != 1 || npar[2] != 1 || npar[3]  != 1) {
        throw GException::model_invalid_parnames(G_READ, xml,
              "LogParabola requires \"Prefactor\" or \"norm\","
              " \"Index\" or \"alpha\", \"Curvature\" or \"beta\""
              " and \"Scale\" or \"Eb\" parameters.");
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element into which model information is written.
 *
 * @exception GException::model_invalid_spectral
 *            Existing XML element is not of type "LogParabola"
 * @exception GException::model_invalid_parnum
 *            Invalid number of model parameters or nodes found in XML element.
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter names found in XML element.
 *
 * Write the LogParabola model information into an XML element. The format
 * of the XML elements is
 *
 *     <spectrum type="LogParabola">
 *       <parameter name="Prefactor" scale=".." value=".." min=".." max=".." free=".."/>
 *       <parameter name="Index"     scale=".." value=".." min=".." max=".." free=".."/>
 *       <parameter name="Curvature" scale=".." value=".." min=".." max=".." free=".."/>
 *       <parameter name="Scale"     scale=".." value=".." min=".." max=".." free=".."/>
 *     </spectrum>
 ***************************************************************************/
void GModelSpectralLogParabola::write(GXmlElement& xml) const
{
    // Set model type
    if (xml.attribute("type") == "") {
        xml.attribute("type", type());
    }

    // Verify model type
    if (xml.attribute("type") != type()) {
        throw GException::model_invalid_spectral(G_WRITE, xml.attribute("type"),
              "Spectral model is not of type \""+type()+"\".");
    }

    // If XML element has 0 nodes then append 4 parameter nodes
    if (xml.elements() == 0) {
        xml.append(GXmlElement("parameter name=\"Prefactor\""));
        xml.append(GXmlElement("parameter name=\"Index\""));
        xml.append(GXmlElement("parameter name=\"Curvature\""));
        xml.append(GXmlElement("parameter name=\"Scale\""));
    }

    // Verify that XML element has exactly 4 parameters
    if (xml.elements() != 4 || xml.elements("parameter") != 4) {
        throw GException::model_invalid_parnum(G_WRITE, xml,
              "LogParabola law model requires exactly 4 parameters.");
    }

    // Set or update model parameter attributes
    int npar[] = {0, 0, 0,0};
    for (int i = 0; i < 4; ++i) {

        // Get parameter element
        GXmlElement* par = xml.element("parameter", i);

        // Handle prefactor
        if (par->attribute("name") == "Prefactor"){
            npar[0]++;
            m_norm.write(*par);
        }

        // Handle index
        else if (par->attribute("name") == "Index") {
            npar[1]++;
            m_index.write(*par);
        }

        // Handle index
        else if (par->attribute("name") == "Curvature"){
              npar[2]++;
              m_curvature.write(*par);
        }

        // Handle pivot energy
        else if (par->attribute("name") == "Scale") {
            m_pivot.write(*par);
            npar[3]++;
        }

    } // endfor: looped over all parameters

    // Check of all required parameters are present
    if (npar[0] != 1 || npar[1] != 1 || npar[2] != 1 || npar[3] != 1) {
        throw GException::model_invalid_parnames(G_WRITE, xml,
              "LogParabola requires \"Prefactor\", \"Index\", \"Curvature\""
              " and \"Scale\" parameters.");
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print LogParabola information
 ***************************************************************************/
std::string GModelSpectralLogParabola::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GModelSpectralLogParabola ===\n");

    // Append information
    result.append(parformat("Number of parameters")+str(size()));
    for (int i = 0; i < size(); ++i) {
        result.append("\n"+m_pars[i]->print());
    }

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GModelSpectralLogParabola::init_members(void)
{
    // Initialise powerlaw normalisation
    m_norm.clear();
    m_norm.name("Prefactor");
    m_norm.unit("ph/cm2/s/MeV");
    m_norm.Scale(1.0);
    m_norm.Value(1.0);          // default: 1.0
    m_norm.Min(0.1);            // min:     0.0
    m_norm.free();
    m_norm.Gradient(0.0);
    m_norm.hasgrad(true);

    // Initialise powerlaw index
    m_index.clear();
    m_index.name("Index");
    m_index.Scale(1.0);
    m_index.Value(-2.0);        // default: -2.0
    m_index.Range(-10.0,+10.0); // range:   [-10,+10]
    m_index.free();
    m_index.Gradient(0.0);
    m_index.hasgrad(true);

    // Initialise curvature
    m_curvature.clear();
    m_curvature.name("Curvature");
    m_curvature.Scale(1.0);
    m_curvature.Value(-0.1);        // default: -2.0
    m_curvature.Range(-10.0,+10.0); // range:   [-10,+10]
    m_curvature.free();
    m_curvature.Gradient(0.0);
    m_curvature.hasgrad(true);

    // Initialise pivot energy
    m_pivot.clear();
    m_pivot.name("PivotEnergy");
    m_pivot.unit("MeV");
    m_pivot.Scale(1.0);
    m_pivot.Value(100.0);       // default: 100
    m_pivot.fix();
    m_pivot.Gradient(0.0);
    m_pivot.hasgrad(true);

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_norm);
    m_pars.push_back(&m_index);
    m_pars.push_back(&m_curvature);
    m_pars.push_back(&m_pivot);

    // Initialise MC cache
    m_mc_emin       = 0.0;
    m_mc_emax       = 0.0;
    m_mc_exponent   = 0.0;
    m_mc_pow_emin   = 0.0;
    m_mc_pow_ewidth = 0.0;
    m_mc_norm       = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model GModelSpectralLogParabola members which should be copied.
 ***************************************************************************/
void GModelSpectralLogParabola::copy_members(const GModelSpectralLogParabola& model)
{
    // Copy members
    m_norm  = model.m_norm;
    m_index = model.m_index;
    m_curvature = model.m_curvature;
    m_pivot = model.m_pivot;

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_norm);
    m_pars.push_back(&m_index);
    m_pars.push_back(&m_curvature);
    m_pars.push_back(&m_pivot);

    // Copy MC cache
    m_mc_emin       = model.m_mc_emin;
    m_mc_emax       = model.m_mc_emax;
    m_mc_exponent   = model.m_mc_exponent;
    m_mc_pow_emin   = model.m_mc_pow_emin;
    m_mc_pow_ewidth = model.m_mc_pow_ewidth;
    m_mc_norm       = model.m_mc_norm;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpectralLogParabola::free_members(void)
{
    // Return
    return;
}

/***********************************************************************//**
 * @brief Update Monte Carlo pre computation cache
 *
 * @param[in] emin Minimum photon energy.
 * @param[in] emax Maximum photon energy.
 *
 * Updates the precomputation cache for Monte Carlo simulations.
 ***************************************************************************/
void GModelSpectralLogParabola::update_mc_cache(const GEnergy& emin,
                                            const GEnergy& emax) const

{
	// Only update if boundaries have changed
	if(emin.MeV() != m_mc_emin || emax.MeV() != m_mc_emax){

        // Store energy boundaries
		m_mc_emin = emin.MeV();
		m_mc_emax = emax.MeV();

		// Find a corresponding power law with the criterion
        // Plaw > LogParabola in the given interval
	    double index_pl;

		// Checking the sign of spectrum curvature
		if (curvature() < 0) {

			// Use the spectral index at the pivot energy of the LogParabola
			index_pl  = index();
			m_mc_norm = norm();
		}
		else {
			// Use a power law which connects the ends of the convex,
            // curved model

			// Plaw index defined by the slope of a straight line in the
            // log-log-plane
			index_pl = std::log(eval(emin)/eval(emax))/std::log(emin.MeV()/emax.MeV());

			// Plaw norm defined such that Plaw = LogParabola at emin
			m_mc_norm = eval(emin)/std::pow(emin.MeV()/pivot(),index_pl);

		}

        // Set precomputation cache
		if (index_pl != -1.0) {
			m_mc_exponent   = index_pl + 1.0;
			m_mc_pow_emin   = std::pow(emin.MeV(), m_mc_exponent);
			m_mc_pow_ewidth = std::pow(emax.MeV(), m_mc_exponent) - m_mc_pow_emin;
		}
		else {
			m_mc_exponent   = 0.0;
			m_mc_pow_emin   = std::log(emin.MeV());
			m_mc_pow_ewidth = std::log(emax.MeV()) - m_mc_pow_emin;
		}

	}

    // Return
    return;
}
