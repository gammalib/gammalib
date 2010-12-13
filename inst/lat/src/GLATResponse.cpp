/***************************************************************************
 *               GLATResponse.cpp  -  GLAST LAT Response class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2010 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GLATResponse.cpp
 * @brief GLATResponse class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iostream>
//#include "GException.hpp"
#include "GLATException.hpp"
#include "GVector.hpp"
#include "GLATResponse.hpp"
#include "GLATPointing.hpp"
#include "GLATEventBin.hpp"
#include "GLATEventCube.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_LOAD                             "load(std::string&, std::string&)"
#define G_GET_FITS_VECTOR        "GLATResponse::get_fits_vector(GFitsTable*,"\
                                                        " std::string&, int)"
#define G_DIFFRSP_BIN     "GLATResponse::diffrsp_bin(GLATEventBin&, GModel&,"\
                                             " GEnergy&, GTime&, GPointing&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Constants __________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                       Constructors/destructors                          =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GLATResponse::GLATResponse(void) : GResponse()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] rsp Response to be copied
 ***************************************************************************/
GLATResponse::GLATResponse(const GLATResponse& rsp) : GResponse(rsp)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(rsp);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GLATResponse::~GLATResponse(void)
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
 * @param[in] rsp Response to be assigned
 ***************************************************************************/
GLATResponse& GLATResponse::operator= (const GLATResponse& rsp)
{
    // Execute only if object is not identical
    if (this != &rsp) {

        // Copy base class members
        this->GResponse::operator=(rsp);

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(rsp);

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
 * @brief Return value of diffuse instrument response function
 *        (units: )
 *
 * @param[in] event Observed event.
 * @param[in] model Source model.
 * @param[in] srcEng True energy of photon.
 * @param[in] srcTime True photon arrival time.
 * @param[in] pnt Instrument pointing information.
 *
 * This method returns the diffuse instrument response for a given event and
 * source model. It handles both event atoms and event bins.
 ***************************************************************************/
double GLATResponse::diffrsp(const GEvent& event, const GModel& model,
                             const GEnergy& srcEng, const GTime& srcTime,
                             const GPointing& pnt) const
{
    // Get IRF value
    double irf = (event.isatom())
                 ? diffrsp_atom((GLATEventAtom&)event, model, srcEng, srcTime, pnt)
                 : diffrsp_bin((GLATEventBin&)event, model, srcEng, srcTime, pnt);

    // Return IRF value
    return irf;
}


/***********************************************************************//**
 * @brief Return livetime fraction
 *        (units: s s-1)
 *
 * @param[in] srcDir True photon direction.
 * @param[in] srcEng True energy of photon.
 * @param[in] srcTime True photon arrival time.
 * @param[in] pnt Pointer to instrument pointing information.
 *
 * @todo Dummy livetime fraction of 0.8.
 ***************************************************************************/
double GLATResponse::live(const GSkyDir& srcDir, const GEnergy& srcEng,
                          const GTime& srcTime, const GPointing& pnt) const
{
    // Dummy
    double live = 0.8;
    
    // Return effective area
    return live;
}


/***********************************************************************//**
 * @brief Return time dispersion
 *        (units: s-1)
 *
 * @param[in] obsTime Observed photon arrival time.
 * @param[in] srcDir True photon direction.
 * @param[in] srcEng True energy of photon.
 * @param[in] srcTime True photon arrival time.
 * @param[in] pnt Pointer to instrument pointing information.
 *
 * The actual implementation of this method assumes no time dispersion,
 * which is equivalent of having a Dirac type time dispersion.
 ***************************************************************************/
double GLATResponse::tdisp(const GTime& obsTime,
                           const GSkyDir& srcDir, const GEnergy& srcEng,
                           const GTime& srcTime, const GPointing& pnt) const
{
    // Dirac time dispersion
    double tdisp = (obsTime == srcTime) ? 1.0 : 0.0;

    // Return time dispersion
    return tdisp;
}


/***********************************************************************//**
 * @brief Return integral over PSF
 *
 * @param[in] srcDir True photon direction.
 * @param[in] srcEng True energy of photon.
 * @param[in] srcTime True photon arrival time.
 * @param[in] pnt Pointer to instrument pointing information.
 * @param[in] roi Region of interest of data selection.
 *
 * @todo Implement integration over ROI.
 ***************************************************************************/
double GLATResponse::npsf(const GSkyDir& srcDir, const GEnergy& srcEng,
                          const GTime& srcTime, const GPointing& pnt,
                          const GRoi& roi) const
{
    // Dummy
    double npsf = 1.0;

    // Return integral
    return npsf;
}


/***********************************************************************//**
 * @brief Return integral over energy dispersion
 *
 * @param[in] srcDir True photon direction.
 * @param[in] srcEng True energy of photon.
 * @param[in] srcTime True photon arrival time.
 * @param[in] pnt Pointer to instrument pointing information.
 * @param[in] ebds Energy boundaries of data selection.
 *
 * @todo Implement integration over energy range.
 ***************************************************************************/
double GLATResponse::nedisp(const GSkyDir& srcDir, const GEnergy& srcEng,
                            const GTime& srcTime, const GPointing& pnt,
                            const GEbounds& ebds) const
{
    // Dummy
    double nedisp = 1.0;

    // Return integral
    return nedisp;
}


/***********************************************************************//**
 * @brief Return integral over time dispersion
 *
 * @param[in] srcDir True photon direction.
 * @param[in] srcEng True energy of photon.
 * @param[in] srcTime True photon arrival time.
 * @param[in] pnt Pointer to instrument pointing information.
 * @param[in] gti Good Time Intervals of data selection.
 *
 * @todo Implement integration over GTIs.
 ***************************************************************************/
double GLATResponse::ntdisp(const GSkyDir& srcDir, const GEnergy& srcEng,
                            const GTime& srcTime, const GPointing& pnt,
                            const GGti& gti) const
{
    // Dummy
    double ntdisp = 1.0;

    // Return integral
    return ntdisp;
}


/***********************************************************************//**
 * @brief Clear instance
***************************************************************************/
void GLATResponse::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GResponse::free_members();

    // Initialise members
    this->GResponse::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
***************************************************************************/
GLATResponse* GLATResponse::clone(void) const
{
    return new GLATResponse(*this);
}


/***********************************************************************//**
 * @brief Load LAT response function
 *
 * @param[in] rspname Name of response (name::front,name::back)
 *
 * @exception GException::rsp_invalid_type
 *            Invalid response type detected.
 *
 * Loads the specified LAT response from the calibration database and
 * performs some pre-calculations for faster response determination.
 *
 * @todo Implement LAT specific exceptions.
 *
 * @todo It seems that the split() method has a bug in that it returns an
 *       empty string in place of the separator.
 ***************************************************************************/
void GLATResponse::load(const std::string& rspname)
{
    // Separate response type
    std::vector<std::string> array = split(rspname, "::");

    // Set response name and type
    if (array.size() == 3) {
        m_rspname = array[0];
        m_rsptype = tolower(array[2]);
    }
    else if (array.size() == 1) {
        m_rspname = array[0];
        m_rsptype = "front";
    }
    else
        throw GException::rsp_invalid_type(G_LOAD, m_rsptype);

    // Convert response type string to section
    if (m_rsptype == "front")
        m_section = 0;
    else if (m_rsptype == "back")
        m_section = 1;
    else
        throw GException::rsp_invalid_type(G_LOAD, m_rsptype);

    // Initialise effective area
    aeff_init();

    // Initialise PSF
    psf_init();

    // Initialise energy dispersion
    edisp_init();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save LAT response
 *
 * @param[in] rspname Response file name.
 ***************************************************************************/
void GLATResponse::save(const std::string& rspname) const
{
    // Allocate FITS file
    GFits file;

    // Open FITS file
    file.open(rspname);

    // Save effective area
    aeff_append(file);

    // Append PSF HDUs
    psf_append(file);

    // Save energy dispersions
    edisp_append(file);

    // Save FITS file
    file.save();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print LAT response information
 ***************************************************************************/
std::string GLATResponse::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GLATResponse ===\n");
    result.append(parformat("Calibration database")+m_caldb+"\n");
    result.append(parformat("Response name")+m_rspname);

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
void GLATResponse::init_members(void)
{
    // Initialise effective area members
    aeff_init_members();

    // Initialise PSF members
    psf_init_members();

    // Initialize energy dispersion members
    edisp_init_members();

    // By default use HANDOFF response database.
    char* handoff = getenv("HANDOFF_IRF_DIR");
    if (handoff != NULL)
        m_caldb.assign(handoff);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] rsp Response.
 ***************************************************************************/
void GLATResponse::copy_members(const GLATResponse& rsp)
{
    // Copy effective area members
    aeff_copy_members(rsp);

    // Copy PSF members
    psf_copy_members(rsp);

    // Copy energy dispersion members
    edisp_copy_members(rsp);

    // Copy remaining members
    m_caldb = rsp.m_caldb;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GLATResponse::free_members(void)
{
    // Free effective area members
    aeff_free_members();

    // Free PSF  members
    psf_free_members();

    // Free energy dispersion members
    edisp_free_members();

    // Free remaining members

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load floating point data from vector column into vector
 *
 * @param[in] hdu FITS table pointer.
 * @param[in] colname Table column name.
 * @param[in] row Table row (staring from 0).
 *
 * Loads the content of a FITS table vector column element into a vector.
 ***************************************************************************/
GVector GLATResponse::get_fits_vector(const GFitsTable* hdu, 
                                      const std::string& colname, 
                                      int row)
{
    // Get pointer to
    GFitsTableCol* ptr = ((GFitsTable*)hdu)->column(colname);
    if (ptr == NULL)
        throw GException::fits_column_not_found(G_GET_FITS_VECTOR, colname);

    // Determine number of entries
    int num = ptr->number();

    // Load data into vector
    GVector data(num);
    for (int i = 0; i < num; ++i)
        data(i) = ptr->real(row,i);

    // Return vector
    return data;
}


/***********************************************************************//**
 * @brief Return value of diffuse instrument response function for event atom
 *
 * @param[in] event Observed event atom.
 * @param[in] model Source model.
 * @param[in] srcEng True energy of photon.
 * @param[in] srcTime True photon arrival time.
 * @param[in] pnt Instrument pointing information.
 ***************************************************************************/
double GLATResponse::diffrsp_atom(const GLATEventAtom& event, const GModel& model,
                                  const GEnergy& srcEng, const GTime& srcTime,
                                  const GPointing& pnt) const
{
    // Initialise IRF with "no response"
    double irf = 0.0;

    // Determine index for diffuse response
    //int inx = event.diffinx(model.name());

    // Return IRF value
    return irf;
}


/***********************************************************************//**
 * @brief Return value of diffuse instrument response function for event bin
 *        
 * @param[in] event Observed event bin.
 * @param[in] model Source model.
 * @param[in] srcEng True energy of photon.
 * @param[in] srcTime True photon arrival time.
 * @param[in] pnt Instrument pointing information.
 *
 * @exception GLATException::diffuse_not_found
 *            Diffuse model not found.
 ***************************************************************************/
double GLATResponse::diffrsp_bin(const GLATEventBin& event, const GModel& model,
                                 const GEnergy& srcEng, const GTime& srcTime,
                                 const GPointing& pnt) const
{
    // Get pointer to event cube
    GLATEventCube* cube = event.cube();

    // Search for diffuse response in event cube
    int idiff = -1;
    for (int i = 0; i < cube->ndiffrsp(); ++i) {
        if (cube->diffname(i) == model.name()) {
            idiff = i;
            break;
        }
    }

    // If diffuse response has not been found then throw an exception
    if (idiff == -1)
        throw GLATException::diffuse_not_found(G_DIFFRSP_BIN, model.name());

    // Get srcmap indices and weighting factors
    GNodeArray* nodes = cube->enodes();
    nodes->set_value(log10(srcEng.MeV()));

    // Compute diffuse response
    GSkymap* map    = cube->diffrsp(idiff);
    double*  pixels = map->pixels() + event.ipix();
    double   irf    = nodes->wgt_left()  * pixels[nodes->inx_left()  * map->npix()] +
                      nodes->wgt_right() * pixels[nodes->inx_right() * map->npix()];

    // Divide by solid angle and ontime since source maps are given in units of
    // counts/pixel/MeV.
    irf /= (event.omega() * event.ontime());

    // Return IRF value
    return irf;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/
