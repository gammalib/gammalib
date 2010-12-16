/***************************************************************************
 *               GLATResponse.cpp  -  Fermi LAT Response class             *
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
#include "GLATResponse.hpp"
#include "GLATPointing.hpp"
#include "GLATEventBin.hpp"
#include "GLATEventCube.hpp"
#include "GLATException.hpp"
#include "GException.hpp"
#include "GVector.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_LOAD                                           "load(std::string&)"
#define G_AEFF                                     "GLATResponse::aeff(int&)"
#define G_PSF                                       "GLATResponse::psf(int&)"
#define G_EDISP                                   "GLATResponse::edisp(int&)"
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
 * @brief Return value of point source instrument response function.
 *
 * @param[in] obsDir Observed photon direction.
 * @param[in] obsEng Observed energy of photon.
 * @param[in] obsTime Observed photon arrival time.
 * @param[in] srcDir True photon direction.
 * @param[in] srcEng True energy of photon.
 * @param[in] srcTime True photon arrival time.
 * @param[in] pnt LAT pointing information.
 *
 * This method implements the default and complete instrument response
 * function (IRF).
 ***************************************************************************/
double GLATResponse::irf(const GInstDir& obsDir, const GEnergy& obsEng,
                         const GTime& obsTime,
                         const GSkyDir& srcDir, const GEnergy& srcEng,
                         const GTime& srcTime, const GPointing& pnt) const
{
    // Get point source IRF components
    double irf  =  live(srcDir,  srcEng, srcTime, pnt);
    //irf        *=  aeff(srcDir,  srcEng, srcTime, pnt);
    //irf        *=   psf(obsDir,  srcDir, srcEng, srcTime, pnt);
    //irf        *= edisp(obsEng,  srcDir, srcEng, srcTime, pnt);
    //irf        *= tdisp(obsTime, srcDir, srcEng, srcTime, pnt);

    // Return IRF value
    return irf;
}


/***********************************************************************//**
 * @brief Return integral of instrument response function.
 *
 * @param[in] srcDir True photon direction.
 * @param[in] srcEng True energy of photon.
 * @param[in] srcTime True photon arrival time.
 * @param[in] pnt Instrument pointing information.
 * @param[in] roi Region of interest of data selection.
 * @param[in] ebds Energy boundaries of data selection.
 * @param[in] gti Good Time Intervals of data selection.
 *
 * This method implements the default and complete integral of the instrument
 * response function (IRF). It may be overwritted by a specific method in the
 * derived class that drops response terms that are not used.
 ***************************************************************************/
double GLATResponse::nirf(const GSkyDir&  srcDir, const GEnergy& srcEng,
                          const GTime& srcTime,  const GPointing& pnt,
                          const GRoi& roi, const GEbounds& ebds,
                          const GGti& gti) const
{
    // Get IRF components
    double nirf  =   live(srcDir, srcEng, srcTime, pnt);
    //nirf        *=   aeff(srcDir, srcEng, srcTime, pnt);
    nirf        *=   npsf(srcDir, srcEng, srcTime, pnt, roi);
    nirf        *= nedisp(srcDir, srcEng, srcTime, pnt, ebds);
    //nirf        *= ntdisp(srcDir, srcEng, srcTime, pnt, gti);

    // Return integrated IRF value
    return nirf;
}


/***********************************************************************//**
 * @brief Return value of diffuse instrument response function
 *        (units: )
 *
 * @param[in] event Observed event.
 * @param[in] model Source model.
 * @param[in] srcEng True energy of photon.
 * @param[in] srcTime True photon arrival time.
 * @param[in] pnt LAT pointing information.
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
 * @brief Load Fermi LAT response from calibration database
 *
 * @param[in] rspname Response name (name/name::front/name::back).
 *
 * @exception GException::rsp_invalid_type
 *            Invalid response type encountered.
 *
 * Loads the specified Fermi LAT response from the calibration database.
 ***************************************************************************/
void GLATResponse::load(const std::string& rspname)
{
    // Free members (cannot call clear() here since this will overwrite the
    // calibration database path m_caldb)
    free_members();

    // Determine response types to be loaded
    std::vector<std::string> array = split(rspname, "::");
    if (array.size() == 1) {
        m_rspname  = array[0];
        m_hasfront = true;
        m_hasback  = true;
    }
    else if (array.size() == 2) {
        m_rspname = array[0];
        if (strip_whitespace(tolower(array[1])) == "front")
            m_hasfront = true;
        else if (strip_whitespace(tolower(array[1])) == "back")
            m_hasback  = true;
        else
            throw GException::rsp_invalid_type(G_LOAD, array[1]);
    }
    else
        throw GException::rsp_invalid_type(G_LOAD, m_rspname);

    // Load front IRF if requested
    if (m_hasfront) {
        std::string aeffname  = m_caldb + "/aeff_"  + m_rspname + "_front.fits";
        std::string psfname   = m_caldb + "/psf_"   + m_rspname + "_front.fits";
        std::string edispname = m_caldb + "/edisp_" + m_rspname + "_front.fits";
        GLATAeff*  aeff  = new GLATAeff(aeffname);
        GLATPsf*   psf   = new GLATPsf(psfname);
        GLATEdisp* edisp = new GLATEdisp(edispname);
        m_aeff.push_back(aeff);
        m_psf.push_back(psf);
        m_edisp.push_back(edisp);
    }

    // Load back IRF if requested
    if (m_hasfront) {
        std::string aeffname  = m_caldb + "/aeff_"  + m_rspname + "_back.fits";
        std::string psfname   = m_caldb + "/psf_"   + m_rspname + "_back.fits";
        std::string edispname = m_caldb + "/edisp_" + m_rspname + "_back.fits";
        GLATAeff*  aeff  = new GLATAeff(aeffname);
        GLATPsf*   psf   = new GLATPsf(psfname);
        GLATEdisp* edisp = new GLATEdisp(edispname);
        m_aeff.push_back(aeff);
        m_psf.push_back(psf);
        m_edisp.push_back(edisp);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save Fermi LAT response in calibration database
 *
 * @param[in] rspname Response name (name/name::front/name::back).
 *
 * Saves Fermi LAT response in a set of FITS files with names
 * aeff_<rspname>_[front/back].fits,
 * psf_<rspname>_[front/back].fits, and
 * edisp_<rspname>_[front/back].fits.
 ***************************************************************************/
void GLATResponse::save(const std::string& rspname) const
{
    // Allocate FITS file
    GFits file;

    // Open FITS file
    file.open(rspname);

    // Save effective area
    //aeff_append(file);

    // Append PSF HDUs
    //psf_append(file);

    // Save energy dispersions
    //edisp_append(file);

    // Save FITS file
    file.save();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return pointer on effective area
 *
 * @param[in] index Response index (starting from 0)
 ***************************************************************************/
GLATAeff* GLATResponse::aeff(const int& index) const
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= m_aeff.size())
        throw GException::out_of_range(G_AEFF, index, 0, m_aeff.size()-1);
    #endif

    // Return pointer on effective area
    return m_aeff[index];
}


/***********************************************************************//**
 * @brief Return pointer on point spread function
 *
 * @param[in] index Response index (starting from 0)
 ***************************************************************************/
GLATPsf* GLATResponse::psf(const int& index) const
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= m_psf.size())
        throw GException::out_of_range(G_PSF, index, 0, m_psf.size()-1);
    #endif

    // Return pointer on point spread function
    return m_psf[index];
}


/***********************************************************************//**
 * @brief Return pointer on energy dispersion
 *
 * @param[in] index Response index (starting from 0)
 ***************************************************************************/
GLATEdisp* GLATResponse::edisp(const int& index) const
{
    // Optionally check if the index is valid
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= m_edisp.size())
        throw GException::out_of_range(G_EDISP, index, 0, m_edisp.size()-1);
    #endif

    // Return pointer on energy dispersion
    return m_edisp[index];
}


/***********************************************************************//**
 * @brief Print LAT response information
 ***************************************************************************/
std::string GLATResponse::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GLATResponse ===");
    result.append("\n"+parformat("Calibration database")+m_caldb);
    result.append("\n"+parformat("Response name")+m_rspname);
    result.append("\n"+parformat("Section(s)"));
    if (m_hasfront && m_hasback)
        result.append("front & back");
    else if (m_hasfront)
        result.append("front");
    else if (m_hasback)
        result.append("back");
    else
        result.append("unknown");
    for (int i = 0; i < size(); ++i) {
        result.append("\n"+m_aeff[i]->print());
        result.append("\n"+m_psf[i]->print());
        result.append("\n"+m_edisp[i]->print());
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
void GLATResponse::init_members(void)
{
    // Initialise members
    m_hasfront = false;
    m_hasback  = false;
    m_aeff.clear();
    m_psf.clear();
    m_edisp.clear();
    
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
    // Copy members
    m_hasfront = rsp.m_hasfront;
    m_hasback  = rsp.m_hasback;
    m_aeff     = rsp.m_aeff;
    m_psf      = rsp.m_psf;
    m_edisp    = rsp.m_edisp;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GLATResponse::free_members(void)
{
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
