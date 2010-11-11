/***************************************************************************
 *        GModelSpatialPtsrc.cpp  -  Spatial point source model class      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2010 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GModelSpatialPtsrc.cpp
 * @brief GModelSpatialPtsrc class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GTools.hpp"
#include "GModelSpatialPtsrc.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_PAR                                  "GModelSpatialPtsrc::par(int)"
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
 * @brief Constructor
 ***************************************************************************/
GModelSpatialPtsrc::GModelSpatialPtsrc(void) : GModelSpatial()
{
    // Initialise private members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor
 *
 * @param[in] dir Position of the point source on the sky.
 ***************************************************************************/
GModelSpatialPtsrc::GModelSpatialPtsrc(const GSkyDir& dir) : GModelSpatial()
{
    // Initialise private members for clean destruction
    init_members();

    // Assign Right Ascension and Declination
    m_ra.value(dir.ra_deg());
    m_dec.value(dir.dec_deg());

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor
 *
 * @param[in] xml XML element containing position information.
 ***************************************************************************/
GModelSpatialPtsrc::GModelSpatialPtsrc(const GXmlElement& xml) : GModelSpatial()
{
    // Initialise private members for clean destruction
    init_members();

    // Read information from XML element
    read(xml);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] model Model from which the instance should be built.
 ***************************************************************************/
GModelSpatialPtsrc::GModelSpatialPtsrc(const GModelSpatialPtsrc& model) :
   GModelSpatial(model)
{
    // Initialise private members for clean destruction
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
 * @param[in] model Model which should be assigned.
 ***************************************************************************/
GModelSpatialPtsrc& GModelSpatialPtsrc::operator= (const GModelSpatialPtsrc& model)
{
    // Execute only if object is not identical
    if (this != &model) {

        // Copy base class members
        this->GModelSpatial::operator=(model);

        // Free members
        free_members();

        // Initialise private members for clean destruction
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
 * @brief Returns pointer to a model parameter
 *
 * @param[in] index Parameter index.
 *
 * @exception GException::out_of_range
 *            Parameter index is out of valid range
 ***************************************************************************/
GModelPar* GModelSpatialPtsrc::par(int index) const
{
    // If index is outside boundary then throw an error
    if (index < 0 || index >= m_npars)
        throw GException::out_of_range(G_PAR, index, 0, m_npars-1);

    // Return parameter pointer
    return m_par[index];
}


/***********************************************************************//**
 * @brief Evaluate function
 *
 * @param[in] srcDir True photon arrival direction.
 *
 * Evaluates the spatial part for a point source model. It implements a delta
 * function with respect to the coordinates of the source.
 *
 * @todo Implement delta function (just returns 1.0 for the moment).
 ***************************************************************************/
double GModelSpatialPtsrc::eval(const GSkyDir& srcDir)
{
    // Set value dependent on source distance
    //double value = (srcDir.dist(ra(), dec()) < 1.0e-6) ? 1.0 : 0.0;
    double value = 1.0;

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
 * source and returns 0 gradients.
 *
 * @todo Implement delta function (just returns 1.0 for the moment).
 * @todo Correctly set parameter gradients (actual version sets gradients to
 * 0).
 ***************************************************************************/
double GModelSpatialPtsrc::eval_gradients(const GSkyDir& srcDir)
{
    // Set value dependent on source distance
    //double value = (srcDir.dist(ra(), dec()) < 1.0e-6) ? 1.0 : 0.0;
    double value = 1.0;

    // Set gradients to 0
    m_ra.gradient(0.0);
    m_dec.gradient(0.0);

    // Return value
    return value;
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
 * is required to have 2 parameters named either 'RA' and 'DEC' or 'GLON'
 * and 'GLAT'.
 ***************************************************************************/
void GModelSpatialPtsrc::read(const GXmlElement& xml)
{
    // Verify that XML element has exactly 2 parameters
    if (xml.elements() != 2 || xml.elements("parameter") != 2)
        throw GException::model_invalid_parnum(G_READ, xml,
              "Point source model requires exactly 2 parameters.");

    // Get pointers on both model parameters
    GXmlElement* par1 = (GXmlElement*)xml.element("parameter", 0);
    GXmlElement* par2 = (GXmlElement*)xml.element("parameter", 1);

    // Get sky direction
    GSkyDir dir;
    if (par1->attribute("name") == "RA" && par2->attribute("name") == "DEC")
        dir.radec_deg(todouble(par1->attribute("value")),
                      todouble(par2->attribute("value")));
    else if (par2->attribute("name") == "RA" && par1->attribute("name") == "DEC")
        dir.radec_deg(todouble(par2->attribute("value")),
                      todouble(par1->attribute("value")));
    else if (par1->attribute("name") == "GLON" && par2->attribute("name") == "GLAT")
        dir.lb_deg(todouble(par1->attribute("value")),
                   todouble(par2->attribute("value")));
    else if (par2->attribute("name") == "GLON" && par1->attribute("name") == "GLAT")
        dir.lb_deg(todouble(par2->attribute("value")),
                   todouble(par1->attribute("value")));
    else
        throw GException::model_invalid_parnames(G_READ, xml,
                          "Require either RA/DEC or GLON/GLAT.");

    // Assign sky direction
    m_ra.value(dir.ra_deg());
    m_dec.value(dir.dec_deg());

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
 * @todo The case that an existing spatial XML element with 'GLON' and 'GLAT'
 * as coordinates is not supported.
 ***************************************************************************/
void GModelSpatialPtsrc::write(GXmlElement& xml) const
{
    // Set model type
    if (xml.attribute("type") == "")
        xml.attribute("type", "SkyDirFunction");

    // Verify model type
    if (xml.attribute("type") != "SkyDirFunction")
        throw GException::model_invalid_spatial(G_WRITE, xml.attribute("type"),
              "Spatial model is not of type \"SkyDirFunction\".");

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


/*==========================================================================
 =                                                                         =
 =                    GModelSpatialPtsrc private methods                   =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GModelSpatialPtsrc::init_members(void)
{
    // Initialise parameters
    m_npars  = 2;
    m_par[0] = &m_ra;
    m_par[1] = &m_dec;

    // Initialise Right Ascension
    m_ra = GModelPar();
    m_ra.name("RA");
    m_ra.unit("deg");
    m_ra.fix();

    // Initialise Declination
    m_dec = GModelPar();
    m_dec.name("DEC");
    m_dec.unit("deg");
    m_dec.fix();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] model GModelSpatialPtsrc members which should be copied.
 ***************************************************************************/
void GModelSpatialPtsrc::copy_members(const GModelSpatialPtsrc& model)
{
    // Copy model parameters (we do not need to copy the rest since it is
    // static)
    m_ra  = model.m_ra;
    m_dec = model.m_dec;

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


/***********************************************************************//**
 * @brief Clone class
***************************************************************************/
GModelSpatialPtsrc* GModelSpatialPtsrc::clone(void) const
{
    return new GModelSpatialPtsrc(*this);
}


/*==========================================================================
 =                                                                         =
 =                        GModelSpatialPtsrc friends                       =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Put model in output stream
 *
 * @param[in] os Output stream into which the model will be dumped
 * @param[in] model Model to be dumped
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GModelSpatialPtsrc& model)
{
    // Put observation in stream
    os << "=== GModelSpatialPtsrc ===" << std::endl;
    os << " Number of parameters ......: " << model.m_npars << std::endl;
    for (int i = 0; i < model.m_npars; ++i) {
        if (i > 0)
            os << std::endl;
        os << " Parameter .................: " << *(model.m_par[i]);
    }

    // Return output stream
    return os;
}
