/***************************************************************************
 *   GModelSpatialRadialProfileGauss.cpp - Gaussian radial profile class   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016 by Juergen Knoedlseder                              *
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
 * @file GModelSpatialRadialProfileGauss.cpp
 * @brief Radial Gaussian profile model class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GMath.hpp"
#include "GXmlElement.hpp"
#include "GModelSpatialRadialProfileGauss.hpp"
#include "GModelSpatialRegistry.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */
const GModelSpatialRadialProfileGauss g_radial_disk_seed;
const GModelSpatialRegistry           g_radial_disk_registry(&g_radial_disk_seed);

/* __ Method name definitions ____________________________________________ */
#define G_READ          "GModelSpatialRadialProfileGauss::read(GXmlElement&)"
#define G_WRITE        "GModelSpatialRadialProfileGauss::write(GXmlElement&)"

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
 *
 * Constructs empty radial Gaussian profile
 ***************************************************************************/
GModelSpatialRadialProfileGauss::GModelSpatialRadialProfileGauss(void) :
                                 GModelSpatialRadialProfile()
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
 * Constructs radial Gaussian profile model by extracting information from
 * an XML element. See the read() method for more information about the
 * expected structure of the XML element.
 ***************************************************************************/
GModelSpatialRadialProfileGauss::GModelSpatialRadialProfileGauss(const GXmlElement& xml) :
                                 GModelSpatialRadialProfile()
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
 * @param[in] model Radial Gaussian profile model.
 *
 * Copies radial Gaussian profile model from another radial profile model.
 ***************************************************************************/
GModelSpatialRadialProfileGauss::GModelSpatialRadialProfileGauss(const GModelSpatialRadialProfileGauss& model) :
                                 GModelSpatialRadialProfile(model)
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
 *
 * Destructs radial Gaussian profile model.
 ***************************************************************************/
GModelSpatialRadialProfileGauss::~GModelSpatialRadialProfileGauss(void)
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
 * @param[in] model Radial Gaussian profile model.
 * @return Radial Gaussian profile model.
 *
 * Assigns radial Gaussian profile model.
 ***************************************************************************/
GModelSpatialRadialProfileGauss& GModelSpatialRadialProfileGauss::operator=(const GModelSpatialRadialProfileGauss& model)
{
    // Execute only if object is not identical
    if (this != &model) {

        // Copy base class members
        this->GModelSpatialRadialProfile::operator=(model);

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
 * @brief Clear radial Gaussian profile model
 *
 * Clears radial Gaussian profile model.
 ***************************************************************************/
void GModelSpatialRadialProfileGauss::clear(void)
{
    // Free class members
    free_members();
    this->GModelSpatialRadialProfile::free_members();
    this->GModelSpatialRadial::free_members();
    this->GModelSpatial::free_members();

    // Initialise members
    this->GModelSpatial::init_members();
    this->GModelSpatialRadial::init_members();
    this->GModelSpatialRadialProfile::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone radial Gaussian profile model
 *
 * @return Pointer to deep copy of radial Gaussian profile model.
 *
 * Returns a deep copy of the radial Gaussian profile model.
 ***************************************************************************/
GModelSpatialRadialProfileGauss* GModelSpatialRadialProfileGauss::clone(void) const
{
    // Clone radial disk model
    return new GModelSpatialRadialProfileGauss(*this);
}


/***********************************************************************//**
 * @brief Return maximum model radius (in radians)
 *
 * @return Maximum model radius (in radians).
 ***************************************************************************/
double GModelSpatialRadialProfileGauss::theta_max(void) const
{
    // Return value
    return (m_sigma.value() * gammalib::deg2rad * 5.0);
}


/***********************************************************************//**
 * @brief Read model from XML element
 *
 * @param[in] xml XML element.
 *
 * Reads the Gaussian radial profile model information from an XML element.
 * The XML element shall have either the format 
 *
 *     <spatialModel type="GaussianProfile">
 *       <parameter name="RA"    scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
 *       <parameter name="DEC"   scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
 *       <parameter name="Sigma" scale="1.0" value="0.45"    min="0.01" max="10"  free="1"/>
 *     </spatialModel>
 *
 * or
 *
 *     <spatialModel type="GaussianProfile">
 *       <parameter name="GLON"  scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
 *       <parameter name="GLAT"  scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
 *       <parameter name="Sigma" scale="1.0" value="0.45"    min="0.01" max="10"  free="1"/>
 *     </spatialModel>
 ***************************************************************************/
void GModelSpatialRadialProfileGauss::read(const GXmlElement& xml)
{
    // Read Gaussian location
    GModelSpatialRadial::read(xml);

    // Read Sigma parameter
    const GXmlElement* par = gammalib::xml_get_par(G_READ, xml, "Sigma");
    m_sigma.read(*par);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element into which model information is written.
 *
 * Writes the Gaussian radial profile model information into an XML element.
 * The XML element will have the format 
 *
 *     <spatialModel type="GaussianProfile">
 *       <parameter name="RA"    scale="1.0" value="83.6331" min="-360" max="360" free="1"/>
 *       <parameter name="DEC"   scale="1.0" value="22.0145" min="-90"  max="90"  free="1"/>
 *       <parameter name="Sigma" scale="1.0" value="0.45"    min="0.01" max="10"  free="1"/>
 *     </spatialModel>
 ***************************************************************************/
void GModelSpatialRadialProfileGauss::write(GXmlElement& xml) const
{
    // Write Gaussian location
    GModelSpatialRadial::write(xml);

    // Write Sigma parameter
    GXmlElement* par = gammalib::xml_need_par(G_WRITE, xml, "Sigma");
    m_sigma.write(*par);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing model information.
 ***************************************************************************/
std::string GModelSpatialRadialProfileGauss::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GModelSpatialRadialProfileGauss ===");

        // Append parameters
        result.append("\n"+gammalib::parformat("Number of parameters"));
        result.append(gammalib::str(size()));
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
void GModelSpatialRadialProfileGauss::init_members(void)
{
    // Initialise Gaussian sigma
    m_sigma.clear();
    m_sigma.name("Sigma");
    m_sigma.unit("deg");
    m_sigma.value(2.778e-4); // 1 arcsec
    m_sigma.min(2.778e-4);   // 1 arcsec
    m_sigma.free();
    m_sigma.scale(1.0);
    m_sigma.gradient(0.0);
    m_sigma.has_grad(false);  // Radial components never have gradients

    // Set parameter pointer(s)
    m_pars.push_back(&m_sigma);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Radial Gaussian model.
 *
 * Copies class members from another radial profile model.
 ***************************************************************************/
void GModelSpatialRadialProfileGauss::copy_members(const GModelSpatialRadialProfileGauss& model)
{
    // Copy members. We do not have to push back the members on the parameter
    // stack as this should have been done by init_members() that was called
    // before. Otherwise we would have sigma twice on the stack.
    m_sigma = model.m_sigma;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpatialRadialProfileGauss::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Radial profile
 *
 * @param[in] theta Angular distance from Gaussian centre (radians).
 * @return Profile value.
 ***************************************************************************/
double GModelSpatialRadialProfileGauss::profile_value(const double& theta) const
{
    // Compute value
    double sigma_rad = m_sigma.value() * gammalib::deg2rad;
    double sigma2    = sigma_rad * sigma_rad;
    double theta2    = theta   * theta;
    double value     = std::exp(-0.5 * theta2 / sigma2) /
                       (gammalib::twopi * sigma2);

    // Return value
    return value;
}
