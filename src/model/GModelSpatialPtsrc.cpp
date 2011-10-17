/***************************************************************************
 *        GModelSpatialPtsrc.cpp  -  Spatial point source model class      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2011 by Jurgen Knodlseder                           *
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
 * @file GModelSpatialPtsrc.cpp
 * @brief Point source spatial model class implementation
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GModelSpatialPtsrc.hpp"
#include "GModelSpatialRegistry.hpp"

/* __ Constants __________________________________________________________ */
const double tolerance =  0.000027777778;    // angular tolerance is 1 arcsec

/* __ Globals ____________________________________________________________ */
const GModelSpatialPtsrc    g_spatial_ptsrc_seed;
const GModelSpatialRegistry g_spatial_ptsrc_registry(&g_spatial_ptsrc_seed);

/* __ Method name definitions ____________________________________________ */
#define G_READ                       "GModelSpatialPtsrc::read(GXmlElement&)"
#define G_WRITE                     "GModelSpatialPtsrc::write(GXmlElement&)"

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
GModelSpatialPtsrc::GModelSpatialPtsrc(void) : GModelSpatial()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor
 *
 * @param[in] dir Sky direction.
 *
 * Creates instance of a point source spatial model using a sky direction.
 ***************************************************************************/
GModelSpatialPtsrc::GModelSpatialPtsrc(const GSkyDir& dir) : GModelSpatial()
{
    // Initialise members
    init_members();

    // Assign direction
    this->dir(dir);

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 *
 * Creates instance of a point source spatial model by extracting information
 * from an XML element. See GModelSpatialPtsrc::read() for more information
 * about the expected structure of the XML element.
 ***************************************************************************/
GModelSpatialPtsrc::GModelSpatialPtsrc(const GXmlElement& xml) : GModelSpatial()
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
 * @param[in] model Point source spatial model.
 ***************************************************************************/
GModelSpatialPtsrc::GModelSpatialPtsrc(const GModelSpatialPtsrc& model) :
   GModelSpatial(model)
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
GModelSpatialPtsrc::~GModelSpatialPtsrc(void)
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
 * @param[in] model Point source spatial model.
 ***************************************************************************/
GModelSpatialPtsrc& GModelSpatialPtsrc::operator= (const GModelSpatialPtsrc& model)
{
    // Execute only if object is not identical
    if (this != &model) {

        // Copy base class members
        this->GModelSpatial::operator=(model);

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
void GModelSpatialPtsrc::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GModelSpatial::free_members();

    // Initialise members
    this->GModelSpatial::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
***************************************************************************/
GModelSpatialPtsrc* GModelSpatialPtsrc::clone(void) const
{
    return new GModelSpatialPtsrc(*this);
}


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] srcDir True photon arrival direction.
 *
 * Evaluates the spatial part for a point source model. It implements a delta
 * function with respect to the coordinates of the source. For numerical
 * reasons a certain tolerance is accepted (typically 0.1 arcsec, i.e. well
 * below the angular resolution of gamma-ray telescopes).
 ***************************************************************************/
double GModelSpatialPtsrc::eval(const GSkyDir& srcDir) const
{
    // Set value dependent on source distance
    double value = (srcDir.dist_deg(dir()) < tolerance) ? 1.0 : 0.0;

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Evaluate function and gradients
 *
 * @param[in] srcDir True photon arrival direction.
 *
 * Evaluates the spatial part for a point source model and the gradient.
 * It implements a delta function with respect to the coordinates of the
 * source.  For numerical reasons a certain tolerance is accepted (typically
 * 0.1 arcsec, i.e. well below the angular resolution of gamma-ray
 * telescopes).
 *
 * This method does not provide valid parameter gradients.
 ***************************************************************************/
double GModelSpatialPtsrc::eval_gradients(const GSkyDir& srcDir) const
{
    // Set value dependent on source distance
    double value = (srcDir.dist_deg(dir()) < tolerance) ? 1.0 : 0.0;

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Returns MC sky direction
 *
 * @param[in] ran Random number generator.
 *
 * This method always returns the directon of the point source.
 ***************************************************************************/
GSkyDir GModelSpatialPtsrc::mc(GRan& ran) const
{
    // Return point source direction
    return dir();
}


/***********************************************************************//**
 * @brief Read model from XML element
 *
 * @param[in] xml XML element containing point source model information.
 *
 * @exception GException::model_invalid_parnum
 *            Invalid number of model parameters found in XML element.
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter names found in XML element.
 *
 * Read the point source information from an XML element. The XML element
 * is required to have 2 parameters named either "RA" and "DEC" or "GLON"
 * and "GLAT".
 ***************************************************************************/
void GModelSpatialPtsrc::read(const GXmlElement& xml)
{
    // Verify that XML element has exactly 2 parameters
    if (xml.elements() != 2 || xml.elements("parameter") != 2)
        throw GException::model_invalid_parnum(G_READ, xml,
              "Point source model requires exactly 2 parameters.");

    // Extract model parameters
    bool has_glon = false;
    bool has_glat = false;
    int  npar[]   = {0, 0};
    for (int i = 0; i < 2; ++i) {

        // Get parameter element
        GXmlElement* par = (GXmlElement*)xml.element("parameter", i);

        // Handle RA/GLON
        if (par->attribute("name") == "RA") {
            m_ra.read(*par);
            npar[0]++;
        }
        else if (par->attribute("name") == "GLON") {
            m_ra.read(*par);
            npar[0]++;
            has_glon = true;
        }

        // Handle DEC/GLAT
        else if (par->attribute("name") == "DEC") {
            m_dec.read(*par);
            npar[1]++;
        }
        else if (par->attribute("name") == "GLAT") {
            m_dec.read(*par);
            npar[1]++;
            has_glat = true;
        }

    } // endfor: looped over all parameters

    // Check if we have to convert GLON/GLAT into RA/DEC
    if (has_glon && has_glat) {
        GSkyDir dir;
        dir.lb_deg(ra(), dec()),
        m_ra.real_value(dir.ra_deg());
        m_dec.real_value(dir.dec_deg());
    }
    else if (has_glon || has_glat) {
        throw GException::model_invalid_parnames(G_READ, xml,
                          "Require either RA/DEC or GLON/GLAT.");
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write model into XML element
 *
 * @param[in] xml XML element into which model information is written.
 *
 * @exception GException::model_invalid_spatial
 *            Existing XML element is not of type 'SkyDirFunction'
 * @exception GException::model_invalid_parnum
 *            Invalid number of model parameters found in XML element.
 * @exception GException::model_invalid_parnames
 *            Invalid model parameter names found in XML element.
 *
 * Write the point source information into an XML element. The XML element
 * has to be of type 'SkyDirFunction' and will have 2 parameter leafs
 * named 'RA' and 'DEC'.
 *
 * @todo The case that an existing spatial XML element with "GLON" and "GLAT"
 *       as coordinates is not supported.
 ***************************************************************************/
void GModelSpatialPtsrc::write(GXmlElement& xml) const
{
    // Set model type
    if (xml.attribute("type") == "")
        xml.attribute("type", type());

    // Verify model type
    if (xml.attribute("type") != type())
        throw GException::model_invalid_spatial(G_WRITE, xml.attribute("type"),
              "Spatial model is not of type \""+type()+"\".");

    // If XML element has 0 nodes then append 2 parameter nodes
    if (xml.elements() == 0) {
        xml.append(new GXmlElement("parameter name=\"RA\""));
        xml.append(new GXmlElement("parameter name=\"DEC\""));
    }

    // Verify that XML element has exactly 2 parameters
    if (xml.elements() != 2 || xml.elements("parameter") != 2)
        throw GException::model_invalid_parnum(G_WRITE, xml,
              "Point source model requires exactly 2 parameters.");

    // Get pointers on both model parameters
    GXmlElement* par1 = (GXmlElement*)xml.element("parameter", 0);
    GXmlElement* par2 = (GXmlElement*)xml.element("parameter", 1);

    // Set or update sky direction
    if (par1->attribute("name") == "RA" && par2->attribute("name") == "DEC") {
        m_ra.write(*par1);
        m_dec.write(*par2);
    }
    else if (par2->attribute("name") == "RA" && par1->attribute("name") == "DEC") {
        m_ra.write(*par2);
        m_dec.write(*par1);
    }
    else
        throw GException::model_invalid_parnames(G_WRITE, xml,
                          "Require RA and DEC parameters.");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print point source information
 ***************************************************************************/
std::string GModelSpatialPtsrc::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GModelSpatialPtsrc ===\n");
    result.append(parformat("Number of parameters")+str(size()));
    for (int i = 0; i < size(); ++i)
        result.append("\n"+m_pars[i]->print());

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Return position of point source
 ***************************************************************************/
GSkyDir GModelSpatialPtsrc::dir(void) const
{
    // Allocate sky direction
    GSkyDir srcDir;

    // Set sky direction
    srcDir.radec_deg(ra(), dec());

    // Return direction
    return srcDir;
}


/***********************************************************************//**
 * @brief Set position of point source
 ***************************************************************************/
void GModelSpatialPtsrc::dir(const GSkyDir& dir)
{
    // Assign Right Ascension and Declination
    m_ra.real_value(dir.ra_deg());
    m_dec.real_value(dir.dec_deg());

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                            Private methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GModelSpatialPtsrc::init_members(void)
{
    // Initialise Right Ascension
    m_ra.clear();
    m_ra.name("RA");
    m_ra.unit("deg");
    m_ra.fix();
    m_ra.scale(1.0);
    m_ra.gradient(0.0);
    m_ra.hasgrad(false);

    // Initialise Declination
    m_dec.clear();
    m_dec.name("DEC");
    m_dec.unit("deg");
    m_dec.fix();
    m_dec.scale(1.0);
    m_dec.gradient(0.0);
    m_dec.hasgrad(false);

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_ra);
    m_pars.push_back(&m_dec);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model Point source spatial model.
 ***************************************************************************/
void GModelSpatialPtsrc::copy_members(const GModelSpatialPtsrc& model)
{
    // Copy members
    m_ra  = model.m_ra;
    m_dec = model.m_dec;

    // Set parameter pointer(s)
    m_pars.clear();
    m_pars.push_back(&m_ra);
    m_pars.push_back(&m_dec);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpatialPtsrc::free_members(void)
{
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Friends                                  =
 =                                                                         =
 ==========================================================================*/
