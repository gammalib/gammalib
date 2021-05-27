/***************************************************************************
 *                 GCTARoi.cpp - CTA region of interest class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2021 by Juergen Knoedlseder                         *
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
 * @file GCTARoi.cpp
 * @brief CTA region of interest class interface implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GEvent.hpp"
#include "GXmlElement.hpp"
#include "GCTARoi.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_READ                                  "GCTARoi::read(GXmlElement&)"
#define G_WRITE                                "GCTARoi::write(GXmlElement&)"
#define G_RADIUS                                   "GCTARoi::radius(double&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */



/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GCTARoi::GCTARoi(void) : GRoi()
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML element constructor
 *
 * @param[in] xml XML element.
 *
 * Constructs region of interest from an XML element.
 ***************************************************************************/
GCTARoi::GCTARoi(const GXmlElement& xml)
{
    // Initialise members
    init_members();

    // Read region of interest from XML element
    read(xml);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Region of interest constructor
 *
 * @param[in] centre Region of interest centre.
 * @param[in] radius Radius of region of interest (degrees).
 ***************************************************************************/
GCTARoi::GCTARoi(const GCTAInstDir& centre, const double& radius) : GRoi()
{
    // Initialise class members
    init_members();

    // Set members
    this->centre(centre);
    this->radius(radius);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] roi Region of interest.
 ***************************************************************************/
GCTARoi::GCTARoi(const GCTARoi& roi) : GRoi(roi)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(roi);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCTARoi::~GCTARoi(void)
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
 * @param[in] roi Region of interest.
 * @return Region of interest.
 ***************************************************************************/
GCTARoi& GCTARoi::operator= (const GCTARoi& roi)
{
    // Execute only if object is not identical
    if (this != &roi) {

        // Copy base class members
        this->GRoi::operator=(roi);

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(roi);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear instance
 ***************************************************************************/
void GCTARoi::clear(void)
{
    // Free members
    free_members();
    this->GRoi::free_members();

    // Initialise private members
    this->GRoi::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
 *
 * @return Pointer to deep copy of region of interest
 ***************************************************************************/
GCTARoi* GCTARoi::clone(void) const
{
    return new GCTARoi(*this);
}


/***********************************************************************//**
 * @brief Check if region of interest contains an event
 *
 * @return True if region of interest contains event, false otherwise
 ***************************************************************************/
bool GCTARoi::contains(const GEvent& event) const
{
    // Initialise flag to non-containment
    bool contains = false;

    // Get pointer to CTA instrument direction
    const GCTAInstDir* dir = dynamic_cast<const GCTAInstDir*>(&event.dir());

    // If instrument direction is a CTA instrument direction then check on
    // containment
    if (dir != NULL) {
        contains = this->contains(*dir);
    }

    // Return containment flag
    return contains;
}


/***********************************************************************//**
 * @brief Check if region of interest contains an instrument direction
 *
 * @return True if region of interest contains instrument direction, false
 * otherwise
 ***************************************************************************/
bool GCTARoi::contains(const GCTAInstDir& dir) const
{
    // Set containment flag
    bool contains = (m_centre.dir().dist_deg(dir.dir()) <= m_radius);

    // Return containment flag
    return contains;
}


/***********************************************************************//**
 * @brief Read region of interest from XML element
 *
 * @param[in] xml XML element.
 *
 * @exception GException::invalid_value
 *            Invalid XML format encountered.
 *
 * Read region of interest from an XML element. The format of the region
 * of interest is
 *
 *     <parameter name="RegionOfInterest" ra="..." dec="..." rad="..."/>
 *
 * The units are deg.
 ***************************************************************************/
void GCTARoi::read(const GXmlElement& xml)
{
    // Clear region of interest
    clear();

    // Get region of interest parameter
    const GXmlElement* par = gammalib::xml_get_par(G_READ, xml, "RegionOfInterest");

    // Extract position attributes
    if (par->has_attribute("ra") &&
        par->has_attribute("dec") &&
        par->has_attribute("rad")) {
        double ra  = gammalib::todouble(par->attribute("ra"));
        double dec = gammalib::todouble(par->attribute("dec"));
        double rad = gammalib::todouble(par->attribute("rad"));
        GSkyDir centre;
        centre.radec_deg(ra,dec);
        this->centre(GCTAInstDir(centre));
        this->radius(rad);
    }
    else {
        std::string msg = "Attributes \"ra\", \"dec\" and/or \"rad\" not found"
                          " in XML parameter \"RegionOfInterest\"."
                          " Please verify the XML format.";
        throw GException::invalid_value(G_READ, msg);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write region of interest into XML element
 *
 * @param[in] xml XML element.
 *
 * Writes region of interest into an XML element. The format of the region
 * of interest is
 *
 *     <parameter name="RegionOfInterest" ra="..." dec="..." rad="..."/>
 *
 * The units are deg.
 *
 * This method does nothing if the region of interest is not valid (signaled
 * by a radius of 0 deg).
 ***************************************************************************/
void GCTARoi::write(GXmlElement& xml) const
{
    // Continue only if the ROI is valid
    if (m_radius > 0.0) {

        // Get parameter
        GXmlElement* par = gammalib::xml_need_par(G_WRITE, xml, "RegionOfInterest");

        // Write attributes
        par->attribute("ra",  gammalib::str(m_centre.dir().ra_deg()));
        par->attribute("dec", gammalib::str(m_centre.dir().dec_deg()));
        par->attribute("rad", gammalib::str(m_radius));

    } // endif: ROI was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print ROI information
 *
 * @param[in] chatter Chattiness.
 * @return String containing ROI information.
 ***************************************************************************/
std::string GCTARoi::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCTARoi ===");

        // Append information
        result.append("\n"+gammalib::parformat("RoI centre"));
        result.append(m_centre.print());
        result.append("\n"+gammalib::parformat("RoI radius"));
        result.append(gammalib::str(m_radius)+" deg");

    } // endif: chatter was not silent

    // Return result
    return result;
}

/***********************************************************************//**
 * @brief Set radius of region of interest
 *
 * @param[in] radius Region of interest radius (degrees).
 *
 * @exception GException::invalid_argument
 *            Non-positive ROI radius encountered
 *
 * Set the radius of the region of interest.
 ***************************************************************************/
void GCTARoi::radius(const double& radius)
{
    // Throw an exception if argument is not valid
    if (radius <= 0.0) {
        std::string msg = "Invalid RoI radius "+gammalib::str(radius)+
                          " encountered. Please specify a strictly"
                          " positive radius.";
        throw GException::invalid_argument(G_RADIUS, msg);
    }

    // Set radius value
    m_radius = radius;

    // Return
    return;
}

/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GCTARoi::init_members(void)
{
    // Initialise members
    m_centre.clear();
    m_radius = 0.0;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] roi Region of interest from which members should be copied.
 ***************************************************************************/
void GCTARoi::copy_members(const GCTARoi& roi)
{
    // Copy attributes
    m_centre = roi.m_centre;
    m_radius = roi.m_radius;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTARoi::free_members(void)
{
    // Return
    return;
}
