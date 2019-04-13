/***************************************************************************
 *  GCTAModelSpatialGaussSpectrum.cpp - Spatial energy dependent Gaussian  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2019 by Juergen Knoedlseder                              *
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
 * @file GCTAModelSpatialGaussSpectrum.cpp
 * @brief Spatial energy dependent Gaussian interface implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include "GTools.hpp"
#include "GMath.hpp"
#include "GCTAInstDir.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GCTAModelSpatialGaussSpectrum.hpp"
#include "GCTAModelSpatialRegistry.hpp"
#include "GModelSpectralRegistry.hpp"
#include "GModelSpectralConst.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GCTAModelSpatialGaussSpectrum g_cta_spatial_gauss_spec_seed;
const GCTAModelSpatialRegistry      g_cta_spatial_gauss_spec_registry(&g_cta_spatial_gauss_spec_seed);

/* __ Method name definitions ____________________________________________ */
#define G_READ            "GCTAModelSpatialGaussSpectrum::read(GXmlElement&)"
#define G_WRITE          "GCTAModelSpatialGaussSpectrum::write(GXmlElement&)"

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
GCTAModelSpatialGaussSpectrum::GCTAModelSpatialGaussSpectrum(void) :
                               GCTAModelSpatial()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Energy-independent Gaussian constructor
 *
 * @param[in] sigma Gaussian width (degrees\f$^2\f$).
 ***************************************************************************/
GCTAModelSpatialGaussSpectrum::GCTAModelSpatialGaussSpectrum(const double& sigma) :
                               GCTAModelSpatial()
{
    // Initialise members
    init_members();

    // Assign sigma
    this->sigma(sigma);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Gaussian spectrum constructor
 *
 * @param[in] sigma Spectrum defining the energy-deependent Gaussian width
 *                  (degrees\f$^2\f$).
 ***************************************************************************/
GCTAModelSpatialGaussSpectrum::GCTAModelSpatialGaussSpectrum(const GModelSpectral& sigma) :
                               GCTAModelSpatial()
{
    // Initialise members
    init_members();

    // Assign sigma spectrum
    this->sigma(sigma);

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 *
 * Creates instance of a spatial energy-dependent Gaussian model by
 * extracting information from an XML element. See
 * GCTAModelSpatialGaussSpectrum::read() for more information about the
 * expected structure of the XML element.
 ***************************************************************************/
GCTAModelSpatialGaussSpectrum::GCTAModelSpatialGaussSpectrum(const GXmlElement& xml) :
                               GCTAModelSpatial()
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
 * @param[in] model Energy-dependent Gaussian model.
 ***************************************************************************/
GCTAModelSpatialGaussSpectrum::GCTAModelSpatialGaussSpectrum(const GCTAModelSpatialGaussSpectrum& model) :
                               GCTAModelSpatial(model)
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
GCTAModelSpatialGaussSpectrum::~GCTAModelSpatialGaussSpectrum(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                               Operators                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] model Energy-dependent Gaussian model.
 ***************************************************************************/
GCTAModelSpatialGaussSpectrum& GCTAModelSpatialGaussSpectrum::operator=(const GCTAModelSpatialGaussSpectrum& model)
{
    // Execute only if object is not identical
    if (this != &model) {

        // Copy base class members
        this->GCTAModelSpatial::operator=(model);

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
 =                            Public methods                               =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear instance
 ***************************************************************************/
void GCTAModelSpatialGaussSpectrum::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GCTAModelSpatial::free_members();

    // Initialise members
    this->GCTAModelSpatial::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
 ***************************************************************************/
GCTAModelSpatialGaussSpectrum* GCTAModelSpatialGaussSpectrum::clone(void) const
{
    return new GCTAModelSpatialGaussSpectrum(*this);
}


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] dir Event direction.
 * @param[in] energy Event energy.
 * @param[in] time Event time.
 * @param[in] gradients Compute gradients?
 * @return Function value
 *
 * Evaluates the energy-dependent Gaussian model for a given event direction.
 * The Gaussian model is defined by
 *
 * \f[f(\theta,E) = \exp \left(-\frac{1}{2}
 *                       \left( \frac{\theta^2}{\sigma(E)} \right)^2 \right)\f]
 *
 * where
 * \f$\theta\f$ is the offset angle (in degrees), and
 * \f$\sigma(E)\f$ is the energy-dependent Gaussian width (in degrees\f$^2\f$).
 *
 * If the @p gradients flag is true the method will also compute the partial
 * derivatives of the parameters. The partial derivative of the Gaussian width
 * is given by
 *
 * \f[\frac{df}{d\sigma_v} = f(\theta) \frac{\theta^4}{\sigma^3} \sigma_s\f]
 *
 * where
 * \f$\sigma_v\f$ is the value part, 
 * \f$\sigma_s\f$ is the scaling part, and 
 * \f$\sigma = \sigma_v \sigma_s\f$.
 *
 * Note that this method implements a function which is unity for
 * \f$\theta=0\f$.
 ***************************************************************************/
double GCTAModelSpatialGaussSpectrum::eval(const GCTAInstDir& dir,
                                           const GEnergy&     energy,
                                           const GTime&       time,
                                           const bool&        gradients) const
{
    // Initialise value
    double value = 0.0;

    // Continue only if a sigma spectrum exists
    if (m_sigma != NULL) {

        // Compute offset angle in degrees
        double offset = dir.theta() * gammalib::rad2deg;

        // Compute energy dependent sigma value in degrees^2
        double sigma = m_sigma->eval(energy, time, gradients);

        // Compute Gaussian value
        double arg  = offset * offset / sigma;
        double arg2 = arg * arg;
        value = std::exp(-0.5 * arg2);

        // Optionally compute partial derivatives
        if (gradients) {

            // Compute partial derivatives of the sigma parameter
            double g_sigma = value * arg2 / sigma;

            // Loop over model parameters
            for (int i = 0; i < m_sigma->size(); ++i) {

                // Get reference to model parameter
                GModelPar& par = m_sigma->operator[](i);

                // Scale parameter gradient
                par.gradient(par.gradient()*g_sigma);

            } // endfor: loop over model parameters

        } // endif: computed partial derivatives

        // Compile option: Check for NaN/Inf
        #if defined(G_NAN_CHECK)
        if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
            std::cout << "*** ERROR: GCTAModelSpatialGaussSpectrum::eval";
            std::cout << "(offset=" << offset << "): NaN/Inf encountered";
            std::cout << " (value=" << value;
            std::cout << ", energy=" << energy;
            std::cout << ", sigma=" << sigma;
            std::cout << ")" << std::endl;
        }
        #endif

    } // endif: sigma spectrum existed

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Read model from XML element
 *
 * @param[in] xml XML element.
 *
 * Read the energy-dependent Gaussian spatial model information from an XML
 * element. The XML element needs to be of the following format:
 *
 *     <spatialModel type="EnergyDependentGaussian">
 *       <sigma type="...">
 *         ...
 *       </sigma>
 *     </spatialModel>
 *
 * where any spectral model type can be provided for the @a sigma tag.
 ***************************************************************************/
void GCTAModelSpatialGaussSpectrum::read(const GXmlElement& xml)
{
    // Clear model
    clear();

    // Get pointers on sigma model components
    const GXmlElement* sigma = xml.element("sigma", 0);

    // Get sigma model
    GModelSpectralRegistry registry;
    m_sigma = registry.alloc(*sigma);

    // Set pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element.
 *
 * @exception GException::invalid_value
 *            Spatial model is not of valid type.
 *
 * Write the energy-dependent Gaussian spatial model information into an XML
 * element. The XML element will be of the following format:
 *
 *     <spatialModel type="EnergyDependentGaussian">
 *       <sigma type="...">
 *         ...
 *       </sigma>
 *     </spatialModel>
 *
 * where the spectral sigma model will be written into the @a sigma tag.
 ***************************************************************************/
void GCTAModelSpatialGaussSpectrum::write(GXmlElement& xml) const
{
    // Set model type
    if (xml.attribute("type") == "") {
        xml.attribute("type", type());
    }

    // Verify model type
    if (xml.attribute("type") != type()) {
        std::string msg = "Spatial model \""+xml.attribute("type")+
                          "\" is not of type \""+type()+"\".";
        throw GException::invalid_value(G_WRITE, msg);
    }

    // Write sigma spectrum if it exists
    if (m_sigma != NULL) {

        // Create new sigma spectrum node
        xml.append(GXmlElement("sigma"));

        // Get new sigma spectrum node
        GXmlElement* sigma = xml.element("sigma", xml.elements("sigma")-1);

        // Write sigma spectrum
        m_sigma->write(*sigma);

    } // endif: sigma spectrum was not NULL

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set sigma spectrum from value
 *
 * @param[in] sigma Constant sigma value.
 ***************************************************************************/
void GCTAModelSpatialGaussSpectrum::sigma(const double& sigma)
{
    // Clear sigma spectrum
    clear();

    // Set constant spectrum
    m_sigma = new GModelSpectralConst(sigma);

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set sigma spectrum
 *
 * @param[in] sigma Sigma spectrum.
 ***************************************************************************/
void GCTAModelSpatialGaussSpectrum::sigma(const GModelSpectral& sigma)
{
    // Clear sigma spectrum
    clear();

    // Set spectrum
    m_sigma = sigma.clone();

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print point source information
 *
 * @param[in] chatter Chattiness.
 * @return String containing point source information.
 ***************************************************************************/
std::string GCTAModelSpatialGaussSpectrum::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCTAModelSpatialGaussSpectrum ===");

        // Determine number of sigma parameters
        int n_sigma = (m_sigma != NULL) ? m_sigma->size() : 0;

        // Append sigma spectrum type
        result.append("\n"+gammalib::parformat("Sigma spectrum"));
        if (n_sigma > 0) {
            result.append("\""+m_sigma->type()+"\"");
        }

        // Append sigma spectrum parameters
        for (int i = 0; i < n_sigma; ++i) {
            result.append("\n"+(m_sigma)[i].print());
        }

    } // endif: chatter was not silent

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                            Private methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GCTAModelSpatialGaussSpectrum::init_members(void)
{
    // Initialise sigma spectrum pointer
    m_sigma = NULL;

    // Clear parameter pointer(s)
    m_pars.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Energy-dependent Gaussian model.
 ***************************************************************************/
void GCTAModelSpatialGaussSpectrum::copy_members(const GCTAModelSpatialGaussSpectrum& model)
{
    // Clone sigma spectrum
    if (model.m_sigma != NULL) {
        m_sigma = model.m_sigma->clone();
    }

    // Set parameter pointers
    set_pointers();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAModelSpatialGaussSpectrum::free_members(void)
{
    // Delete sigma spectrum
    if (m_sigma != NULL) delete m_sigma;

    // Signal free pointer
    m_sigma = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set pointers
 *
 * Set pointers to all model parameters. The pointers are stored in a vector
 * that is member of the GModelData base class.
 ***************************************************************************/
void GCTAModelSpatialGaussSpectrum::set_pointers(void)
{
    // Clear parameters
    m_pars.clear();

    // Determine number of sigma parameters
    int n_sigma = (m_sigma != NULL) ? m_sigma->size() : 0;

    // Gather sigma parameters if there are some
    if (n_sigma > 0) {
        for (int i = 0; i < n_sigma; ++i) {
            m_pars.push_back(&((*m_sigma)[i]));
        }
    }

    // Return
    return;
}
