/***************************************************************************
 *     GModelSpatialPointSource.cpp - Spatial point source model class     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2013 by Juergen Knoedlseder                         *
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
 * @file GModelSpatialPointSource.cpp
 * @brief Point source spatial model class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GModelSpatialPointSource.hpp"
#include "GModelSpatialRegistry.hpp"

/* __ Constants __________________________________________________________ */
const double tolerance =  0.000027777778;    // angular tolerance is 1 arcsec

/* __ Globals ____________________________________________________________ */
const GModelSpatialPointSource g_spatial_ptsrc_seed;
const GModelSpatialRegistry    g_spatial_ptsrc_registry(&g_spatial_ptsrc_seed);

/* __ Method name definitions ____________________________________________ */
#define G_READ                 "GModelSpatialPointSource::read(GXmlElement&)"
#define G_WRITE               "GModelSpatialPointSource::write(GXmlElement&)"

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
GModelSpatialPointSource::GModelSpatialPointSource(void) : GModelSpatial()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Sky direction constructor
 *
 * @param[in] dir Sky direction.
 *
 * Construct a point source spatial model from a sky direction.
 ***************************************************************************/
GModelSpatialPointSource::GModelSpatialPointSource(const GSkyDir& dir) :
                          GModelSpatial()
{
    // Initialise members
    init_members();

    // Assign direction
    this->dir(dir);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Value constructor
 *
 * @param[in] ra Right Ascencion of model centre.
 * @param[in] dec Declination of model centre.
 *
 * Construct a point source spatial model from the Right Ascension and
 * Declination of the model centre.
 ***************************************************************************/
GModelSpatialPointSource::GModelSpatialPointSource(const double& ra,
                                                   const double& dec) :
                          GModelSpatial()
{
    // Initialise members
    init_members();

    // Set values
    m_ra.value(ra);
    m_dec.value(dec);

    // Return
    return;
}


/***********************************************************************//**
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 *
 * Construct a point source spatial model by extracting information from an
 * XML element. See the read() method for more information about the expected
 * structure of the XML element.
 ***************************************************************************/
GModelSpatialPointSource::GModelSpatialPointSource(const GXmlElement& xml) :
                          GModelSpatial()
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
GModelSpatialPointSource::GModelSpatialPointSource(const GModelSpatialPointSource& model) :
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
GModelSpatialPointSource::~GModelSpatialPointSource(void)
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
 * @return Point source spatial model.
 ***************************************************************************/
GModelSpatialPointSource& GModelSpatialPointSource::operator=(const GModelSpatialPointSource& model)
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
 * @brief Clear point source model
 ***************************************************************************/
void GModelSpatialPointSource::clear(void)
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
 * @brief Clone point source model
 *
 * @return Pointer to deep copy of point source model.
 ***************************************************************************/
GModelSpatialPointSource* GModelSpatialPointSource::clone(void) const
{
    // Clone point source model
    return new GModelSpatialPointSource(*this);
}


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] photon Incident photon.
 * @return Model value.
 *
 * Evaluates the spatial part for a point source model. It implements a delta
 * function with respect to the coordinates of the source. For numerical
 * reasons, a certain tolerance is accepted (typically 0.1 arcsec, i.e. well
 * below the angular resolution of gamma-ray telescopes).
 ***************************************************************************/
double GModelSpatialPointSource::eval(const GPhoton& photon) const
{
    // Set value dependent on source distance
    double value = (photon.dir().dist_deg(dir()) < tolerance) ? 1.0 : 0.0;

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Evaluate function and gradients
 *
 * @param[in] photon Incident photon.
 * @return Model value.
 *
 * Evaluates the spatial part for a point source model and the gradient.
 * It implements a delta function with respect to the coordinates of the
 * source.  For numerical reasons a certain tolerance is accepted (typically
 * 0.1 arcsec, i.e. well below the angular resolution of gamma-ray
 * telescopes).
 *
 * This method does not provide valid parameter gradients.
 ***************************************************************************/
double GModelSpatialPointSource::eval_gradients(const GPhoton& photon) const
{
    // Return value
    return (eval(photon));
}


/***********************************************************************//**
 * @brief Returns MC sky direction
 *
 * @param[in] energy Photon energy.
 * @param[in] time Photon arrival time.
 * @param[in,out] ran Random number generator.
 * @return Sky direction.
 *
 * Draws an arbitrary sky direction for the point source model. As the point
 * source is a point in the sky, the methods always returns the directon of
 * the point source.
 ***************************************************************************/
GSkyDir GModelSpatialPointSource::mc(const GEnergy& energy,
                                     const GTime&   time,
                                     GRan&          ran) const
{
    // Return sky direction
    return (dir());
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
 * Read the point source information from an XML element with the following
 * format
 *
 *     <spatialModel type="SkyDirFunction">
 *       <parameter free="0" max="360" min="-360" name="RA" scale="1" value="83.6331" />
 *       <parameter free="0" max="90" min="-90" name="DEC" scale="1" value="22.0145" />
 *     </spatialModel>
 *
 * or
 *
 *     <spatialModel type="SkyDirFunction">
 *       <parameter free="0" max="360" min="-360" name="GLON" scale="1" value="83.6331" />
 *       <parameter free="0" max="90" min="-90" name="GLAT" scale="1" value="22.0145" />
 *     </spatialModel>
 *
 ***************************************************************************/
void GModelSpatialPointSource::read(const GXmlElement& xml)
{
    // Verify that XML element has exactly 2 parameters
    if (xml.elements() != 2 || xml.elements("parameter") != 2) {
        throw GException::model_invalid_parnum(G_READ, xml,
              "Point source model requires exactly 2 parameters.");
    }

    // Extract model parameters
    bool has_glon = false;
    bool has_glat = false;
    int  npar[]   = {0, 0};
    for (int i = 0; i < 2; ++i) {

        // Get parameter element
        const GXmlElement* par = xml.element("parameter", i);

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
        m_ra.value(dir.ra_deg());
        m_dec.value(dir.dec_deg());
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
 * Write the point source information into an XML element with the following
 * format
 *
 *     <spatialModel type="SkyDirFunction">
 *       <parameter free="0" max="360" min="-360" name="RA" scale="1" value="83.6331" />
 *       <parameter free="0" max="90" min="-90" name="DEC" scale="1" value="22.0145" />
 *     </spatialModel>
 *
 * @todo The case that an existing spatial XML element with "GLON" and "GLAT"
 *       as coordinates is not supported.
 ***************************************************************************/
void GModelSpatialPointSource::write(GXmlElement& xml) const
{
    // Set model type
    if (xml.attribute("type") == "") {
        xml.attribute("type", type());
    }

    // Verify model type
    if (xml.attribute("type") != type()) {
        throw GException::model_invalid_spatial(G_WRITE, xml.attribute("type"),
              "Spatial model is not of type \""+type()+"\".");
    }

    // If XML element has 0 nodes then append 2 parameter nodes
    if (xml.elements() == 0) {
        xml.append(GXmlElement("parameter name=\"RA\""));
        xml.append(GXmlElement("parameter name=\"DEC\""));
    }

    // Verify that XML element has exactly 2 parameters
    if (xml.elements() != 2 || xml.elements("parameter") != 2) {
        throw GException::model_invalid_parnum(G_WRITE, xml,
              "Point source model requires exactly 2 parameters.");
    }

    // Get pointers on both model parameters
    GXmlElement* par1 = xml.element("parameter", 0);
    GXmlElement* par2 = xml.element("parameter", 1);

    // Set or update sky direction
    if (par1->attribute("name") == "RA" && par2->attribute("name") == "DEC") {
        m_ra.write(*par1);
        m_dec.write(*par2);
    }
    else if (par2->attribute("name") == "RA" && par1->attribute("name") == "DEC") {
        m_ra.write(*par2);
        m_dec.write(*par1);
    }
    else {
        throw GException::model_invalid_parnames(G_WRITE, xml,
                          "Require RA and DEC parameters.");
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print point source information
 *
 * @return String containing point source information.
 ***************************************************************************/
std::string GModelSpatialPointSource::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GModelSpatialPointSource ===");

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


/***********************************************************************//**
 * @brief Return position of point source
 *
 * @return Point source sky direction.
 *
 * Returns the sky direction of the point source.
 ***************************************************************************/
GSkyDir GModelSpatialPointSource::dir(void) const
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
 *
 * Sets the sky direction of the point source.
 ***************************************************************************/
void GModelSpatialPointSource::dir(const GSkyDir& dir)
{
    // Assign Right Ascension and Declination
    m_ra.value(dir.ra_deg());
    m_dec.value(dir.dec_deg());

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
void GModelSpatialPointSource::init_members(void)
{
    // Initialise Right Ascension
    m_ra.clear();
    m_ra.name("RA");
    m_ra.unit("deg");
    m_ra.fix();
    m_ra.scale(1.0);
    m_ra.gradient(0.0);
    m_ra.has_grad(false);

    // Initialise Declination
    m_dec.clear();
    m_dec.name("DEC");
    m_dec.unit("deg");
    m_dec.fix();
    m_dec.scale(1.0);
    m_dec.gradient(0.0);
    m_dec.has_grad(false);

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
void GModelSpatialPointSource::copy_members(const GModelSpatialPointSource& model)
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
void GModelSpatialPointSource::free_members(void)
{
    // Return
    return;
}
