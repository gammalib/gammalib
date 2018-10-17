/***************************************************************************
 *      GCTAModelSpatialGradient.cpp - Spatial gradient CTA model class    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2018 by Juergen Knoedlseder                              *
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
 * @file GCTAModelSpatialGradient.cpp
 * @brief Spatial gradient CTA interface definition
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include "GException.hpp"
#include "GTools.hpp"
#include "GMath.hpp"
#include "GRan.hpp"
#include "GCTAModelSpatialGradient.hpp"
#include "GCTAModelSpatialRegistry.hpp"
#include "GCTAInstDir.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GCTAModelSpatialGradient g_cta_spatial_gradient_seed;
const GCTAModelSpatialRegistry g_cta_spatial_gradient_registry(&g_cta_spatial_gradient_seed);

/* __ Method name definitions ____________________________________________ */
#define G_READ                 "GCTAModelSpatialGradient::read(GXmlElement&)"
#define G_WRITE               "GCTAModelSpatialGradient::write(GXmlElement&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
//#define G_DEBUG_MC                                     //!< Debug MC method

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GCTAModelSpatialGradient::GCTAModelSpatialGradient(void) : GCTAModelSpatial()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 *
 * Creates instance of a spatial gradient model by extracting information
 * from an XML element. See GCTAModelSpatialGradient::read() for more
 * information about the expected structure of the XML element.
 ***************************************************************************/
GCTAModelSpatialGradient::GCTAModelSpatialGradient(const GXmlElement& xml) : GCTAModelSpatial()
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
 * @param[in] model Spatial gradient model.
 ***************************************************************************/
GCTAModelSpatialGradient::GCTAModelSpatialGradient(const GCTAModelSpatialGradient& model) :
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
GCTAModelSpatialGradient::~GCTAModelSpatialGradient(void)
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
 * @param[in] model Spatial gradient model.
 ***************************************************************************/
GCTAModelSpatialGradient& GCTAModelSpatialGradient::operator=(const GCTAModelSpatialGradient& model)
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
void GCTAModelSpatialGradient::clear(void)
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
GCTAModelSpatialGradient* GCTAModelSpatialGradient::clone(void) const
{
    return new GCTAModelSpatialGradient(*this);
}


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] dir Event direction.
 * @param[in] energy Event energy (not used).
 * @param[in] time Event time (not used).
 * @param[in] gradients Compute gradients?
 * @return Function value
 *
 * Evaluates the spatial gradient model for a given event direction. The
 * energy and time of the event are not used.
 *
 * The spatial gradient model is defined as
 *
 * \f[f(x,y) = 1 + g_x \times x + g_y \times y\f]
 *
 * where
 * \f$x\f$ is x direction,
 * \f$y\f$ is y direction,
 * \f$g_x\f$ is the spatial gradient in the x direction, and
 * \f$g_y\f$ is the spatial gradient in the y direction.
 *
 * If the @p gradients flag is true the method will also compute the partial
 * derivatives of the parameters. The partial derivative of the spatial
 * gradient model are given by
 *
 * \f[ \frac{df}{dg_{xv}} = g_{xs} x\f]
 *
 * and
 *
 * \f[ \frac{df}{dg_{yv}} = g_{ys} y\f]
 *
 * where
 * \f$g_{xv}\f$ is the value part and \f$g_{xs}\f$ is the scaling part of
 * gradient in x, and
 * \f$g_{yv}\f$ is the value part and \f$g_{ys}\f$ is the scaling part of
 * gradient in y.
 ***************************************************************************/
double GCTAModelSpatialGradient::eval(const GCTAInstDir& dir,
                                      const GEnergy&     energy,
                                      const GTime&       time,
                                      const bool&        gradients) const
{
    // Get detx and dety in degrees
    double detx = dir.detx() * gammalib::rad2deg;
    double dety = dir.dety() * gammalib::rad2deg;

    // Compute value
    double value = 1.0 + m_detx_gradient.value() * detx +
                         m_dety_gradient.value() * dety;

    // Optionally compute partial derivatives
    if (gradients) {
        m_detx_gradient.factor_gradient(m_detx_gradient.scale() * detx);
        m_dety_gradient.factor_gradient(m_dety_gradient.scale() * dety);
    }

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(value) || gammalib::is_infinite(value)) {
        std::cout << "*** ERROR: GCTAModelSpatialGradient::eval";
        std::cout << "(detx=" << detx << ", dety=" << dety;
        std::cout << "): NaN/Inf encountered";
        std::cout << " (value=" << value;
        std::cout << ")" << std::endl;
    }
    #endif

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Returns MC instrument direction
 *
 * @param[in] energy Event energy (not used).
 * @param[in] time Event time (not used).
 * @param[in,out] ran Random number generator.
 * @return CTA instrument direction.
 *
 * @todo Implement method
 ***************************************************************************/
GCTAInstDir GCTAModelSpatialGradient::mc(const GEnergy& energy,
                                         const GTime&   time,
                                         GRan& ran) const
{
    // Set instrument direction
    GCTAInstDir dir;

    // Return instrument direction
    return dir;
}


/***********************************************************************//**
 * @brief Read model from XML element
 *
 * @param[in] xml XML element.
 *
 * Read the gradient spatial model information from an XML element.
 ***************************************************************************/
void GCTAModelSpatialGradient::read(const GXmlElement& xml)
{
    // Get parameter pointers
    const GXmlElement* detx = gammalib::xml_get_par(G_READ, xml, m_detx_gradient.name());
    const GXmlElement* dety = gammalib::xml_get_par(G_READ, xml, m_dety_gradient.name());

    // Read parameters
    m_detx_gradient.read(*detx);
    m_dety_gradient.read(*dety);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element.
 *
 * @exceptino GException::invalid_value
 *            Spatial model is not of valid type.
 *
 * Write the gradient spatial model information into an XML element.
 ***************************************************************************/
void GCTAModelSpatialGradient::write(GXmlElement& xml) const
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

    // Get XML parameters
    GXmlElement* detx = gammalib::xml_need_par(G_WRITE, xml, m_detx_gradient.name());
    GXmlElement* dety = gammalib::xml_need_par(G_WRITE, xml, m_dety_gradient.name());

    // Write parameters
    m_detx_gradient.write(*detx);
    m_dety_gradient.write(*dety);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print point source information
 *
 * @param[in] chatter Chattiness.
 * @return String containing point source information.
 ***************************************************************************/
std::string GCTAModelSpatialGradient::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCTAModelSpatialGradient ===");

        // Append information
        result.append("\n"+gammalib::parformat("Number of parameters") +
                      gammalib::str(size()));
        for (int i = 0; i < size(); ++i) {
            result.append("\n"+m_pars[i]->print(chatter));
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
void GCTAModelSpatialGradient::init_members(void)
{
    // Initialise detx gradient
    m_detx_gradient.clear();
    m_detx_gradient.name("Grad_DETX");
    m_detx_gradient.unit("deg^-1");
    m_detx_gradient.value(0.0);
    m_detx_gradient.free();
    m_detx_gradient.scale(1.0);
    m_detx_gradient.gradient(0.0);
    m_detx_gradient.has_grad(true);

    // Initialise dety gradient
    m_dety_gradient.clear();
    m_dety_gradient.name("Grad_DETY");
    m_dety_gradient.unit("deg^-1");
    m_dety_gradient.value(0.0);
    m_dety_gradient.free();
    m_dety_gradient.scale(1.0);
    m_dety_gradient.gradient(0.0);
    m_dety_gradient.has_grad(true);

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_detx_gradient);
    m_pars.push_back(&m_dety_gradient);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Radial Gaussian model.
 ***************************************************************************/
void GCTAModelSpatialGradient::copy_members(const GCTAModelSpatialGradient& model)
{
    // Copy members
    m_detx_gradient = model.m_detx_gradient;
    m_dety_gradient = model.m_dety_gradient;

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_detx_gradient);
    m_pars.push_back(&m_dety_gradient);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAModelSpatialGradient::free_members(void)
{
    // Return
    return;
}
