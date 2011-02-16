/***************************************************************************
 *               GLATResponse.cpp  -  Fermi LAT Response class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2011 by Jurgen Knodlseder                           *
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
 * @brief Fermi LAT Response class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <unistd.h>           // access() function
#include <cstdlib>            // std::getenv() function
#include <string>
#include "GLATResponse.hpp"
#include "GLATObservation.hpp"
#include "GLATEventAtom.hpp"
#include "GLATEventBin.hpp"
#include "GLATEventCube.hpp"
#include "GLATException.hpp"
#include "GException.hpp"
#include "GModelSpatialPtsrc.hpp"
#include "GFits.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_CALDB                           "GLATResponse::caldb(std::string&)"
#define G_LOAD                             "GLATResponse::load(std::string&)"
#define G_AEFF                                     "GLATResponse::aeff(int&)"
#define G_PSF                                       "GLATResponse::psf(int&)"
#define G_EDISP                                   "GLATResponse::edisp(int&)"
#define G_IRF_ATOM     "GLATResponse::irf(GLATEventAtom&, GModel&, GEnergy&,"\
                                                     "GTime&, GObservation&)"
#define G_IRF_BIN       "GLATResponse::irf(GLATEventBin&, GModel&, GEnergy&,"\
                                                     "GTime&, GObservation&)"
#define G_NPRED              "GLATResponse::npred(GModel&, GEnergy&, GTime&,"\
                                                            " GObservation&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
#define G_DUMP_MEAN_PSF  1                    //!< Dump mean PSF allocation
#define G_DEBUG_MEAN_PSF 0                    //!< Debug mean PSF computation

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
 * @param[in] rsp Response.
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
 * @param[in] rsp Response.
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
 * @brief Set path to the calibration database
 *
 * @param[in] caldb Path to calibration database
 *
 * @exception GException::caldb_not_found
 *            Calibration database repository not found.
 *
 * This default method simply checks if the calibration database directory
 * exists. If the directory exists, the path will be stored. No checking is
 * implemented that checks for the consistency of the calibration database.
 *
 * @todo Implement a GCalDB class that handles any calibration database
 *       issues. GCalDB may be an abstract class for which instrument
 *       specific methods are implement to handle any instrument specific
 *       IRF database issues. 
 ***************************************************************************/
void GLATResponse::caldb(const std::string& caldb)
{
    // Check if calibration database directory is accessible
    if (access(caldb.c_str(), R_OK) != 0)
        throw GException::caldb_not_found(G_CALDB, caldb);
    
    // Store the path to the calibration database
    m_caldb = caldb;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return value of model instrument response function
 *
 * @param[in] event Event.
 * @param[in] model Source model.
 * @param[in] srcEng True energy of photon.
 * @param[in] srcTime True photon arrival time.
 * @param[in] obs Observation.
 *
 * This method returns the response of the instrument to a specific source
 * model. The method handles both event atoms and event bins.
 ***************************************************************************/
double GLATResponse::irf(const GEvent& event, const GModelSky& model,
                         const GEnergy& srcEng, const GTime& srcTime,
                         const GObservation& obs) const
{
    // Get IRF value
    double rsp;
    if (event.isatom())
        rsp = irf(static_cast<const GLATEventAtom&>(event), model,
                  srcEng, srcTime, obs);
    else
        rsp = irf(static_cast<const GLATEventBin&>(event), model,
                  srcEng, srcTime, obs);

    // Return IRF value
    return rsp;
}


/***********************************************************************//**
 * @brief Return value of model IRF for event atom
 *
 * @param[in] event Event atom.
 * @param[in] model Source model.
 * @param[in] srcEng True energy of photon.
 * @param[in] srcTime True photon arrival time.
 * @param[in] obs Observations.
 *
 * @todo Not yet implemented.
 ***************************************************************************/
double GLATResponse::irf(const GLATEventAtom& event, const GModelSky& model,
                         const GEnergy& srcEng, const GTime& srcTime,
                         const GObservation& obs) const
{
    // Initialise IRF with "no response"
    double irf = 0.0;

    // Dump warning that integration is not yet implemented
    throw GException::feature_not_implemented(G_IRF_ATOM);

    // Return IRF value
    return irf;
}


/***********************************************************************//**
 * @brief Return value of model IRF for event bin
 *        
 * @param[in] event Event bin.
 * @param[in] model Source model.
 * @param[in] srcEng True energy of photon.
 * @param[in] srcTime True photon arrival time.
 * @param[in] obs Observation.
 *
 * @exception GLATException::diffuse_not_found
 *            Diffuse model not found.
 *
 * This method first searches for a corresponding source map, and if found,
 * computes the response from the source map. If no source map is present
 * it checks if the source is a point source. If this is the case a mean
 * PSF is allocated for the source and the response is computed from the
 * mean PSF. Otherwise an GLATException::diffuse_not_found exception is
 * thrown.
 *
 * @todo Extract event cube from observation. We do not need the cube
 *       pointer in the event anymore.
 * @todo Implement more efficient method for response search (name search is
 *       not very rapid).
 * @todo Instead of calling "offset = event.dir().dist_deg(srcDir)" we can
 *       precompute and store for each PSF the offsets. This should save
 *       quite some time since the distance computation is time
 *       consuming.
 ***************************************************************************/
double GLATResponse::irf(const GLATEventBin& event, const GModelSky& model,
                         const GEnergy& srcEng, const GTime& srcTime,
                         const GObservation& obs) const
{
    // Initialise response value
    double rsp = 0.0;

    // Get pointer to event cube
    GLATEventCube* cube = event.cube();

    // Get pointer on point source spatial model
    GModelSpatialPtsrc* ptsrc = dynamic_cast<GModelSpatialPtsrc*>(model.spatial());

    // Search for diffuse response in event cube
    int idiff = -1;
    for (int i = 0; i < cube->ndiffrsp(); ++i) {
        if (cube->diffname(i) == model.name()) {
            idiff = i;
            break;
        }
    }

    // If diffuse response has been found then get response from source map
    if (idiff != -1) {

        // Get srcmap indices and weighting factors
        GNodeArray nodes = cube->enodes();
        nodes.set_value(srcEng.log10MeV());

        // Compute diffuse response
        GSkymap* map    = cube->diffrsp(idiff);
        double*  pixels = map->pixels() + event.ipix();
        rsp             = nodes.wgt_left()  * pixels[nodes.inx_left()  * map->npix()] +
                          nodes.wgt_right() * pixels[nodes.inx_right() * map->npix()];

        // Divide by solid angle and ontime since source maps are given in units of
        // counts/pixel/MeV.
        rsp /= (event.omega() * event.ontime());

    } // endelse: diffuse response was present

    // ... otherwise check if model is a point source. If this is true
    // then return response from mean PSF
    if ((idiff == -1 || m_force_mean) && ptsrc != NULL) {

        // Search for mean PSF
        int ipsf = -1;
        for (int i = 0; i < m_ptsrc.size(); ++i) {
            if (m_ptsrc[i]->name() == model.name()) {
                ipsf = i;
                break;
            }
        }

        // If mean PSF has not been found then create it now
        if (ipsf == -1) {

            // Allocate new mean PSF
            GLATMeanPsf* psf = 
                new GLATMeanPsf(ptsrc->dir(), static_cast<const GLATObservation&>(obs));

            // Set source name
            psf->name(model.name());

            // Push mean PSF on stack
            ((GLATResponse*)this)->m_ptsrc.push_back(psf);

            // Set index of mean PSF
            ipsf = m_ptsrc.size()-1;

            // Debug option: dump mean PSF
            #if G_DUMP_MEAN_PSF 
            std::cout << "Added new mean PSF \""+model.name() << "\"" << std::endl;
            std::cout << *psf << std::endl;
            #endif

        } // endif: created new mean PSF

        // Get PSF value
        GSkyDir srcDir   = m_ptsrc[ipsf]->dir();
        double  offset   = event.dir().dist_deg(srcDir);
        double  mean_psf = (*m_ptsrc[ipsf])(offset, srcEng.log10MeV()) / (event.ontime());

        // Debug option: compare mean PSF to diffuse response
        #if G_DEBUG_MEAN_PSF
        std::cout << "Energy=" << srcEng.MeV();
        std::cout << " MeanPsf=" << mean_psf;
        std::cout << " DiffusePsf=" << rsp;
        std::cout << " ratio(Mean/Diffuse)=" << mean_psf/rsp;
        if (mean_psf/rsp < 0.99)
            std::cout << " <<< (1%)";
        else if (mean_psf/rsp > 1.01)
            std::cout << " >>> (1%)";
        std::cout << std::endl;
        #endif

        // Set response
        rsp = mean_psf;

    } // endif: model was point source

    // ... otherwise throw an exception
    if ((idiff == -1) && ptsrc == NULL)
        throw GLATException::diffuse_not_found(G_IRF_BIN, model.name());

    // Return IRF value
    return rsp;
}


/***********************************************************************//**
 * @brief Return integral of instrument response function.
 *
 * @param[in] model Source model.
 * @param[in] srcEng True energy of photon.
 * @param[in] srcTime True photon arrival time.
 * @param[in] obs Observation.
 *
 * @todo Not yet implemented.
 ***************************************************************************/
double GLATResponse::npred(const GModelSky& model, const GEnergy& srcEng,
                           const GTime& srcTime,
                           const GObservation& obs) const
{
    // Initialise
    double npred = 0.0;

    // Notify that method is not yet implemented
    throw GException::feature_not_implemented(G_NPRED);

    // Return integrated IRF value
    return npred;
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
 *
 * @todo Add a more generic calibration database interface. For this we need
 *       first to develop a clear vision about how the calibration database
 *       is organized. Probably it is best to introduce a generic base class
 *       for the calibration database and each instrument implements then
 *       it's specific interface.
 ***************************************************************************/
void GLATResponse::load(const std::string& rspname)
{
    // Save calibration database name
    std::string caldb = m_caldb;

    // Clear instance
    clear();

    // Restore calibration database name
    m_caldb = caldb;

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
    if (m_hasback) {
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
 *
 * @todo Not yet implemented.
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
 * @param[in] index Response index (starting from 0).
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
 * @param[in] index Response index (starting from 0).
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
 * @param[in] index Response index (starting from 0).
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
    for (int i = 0; i < m_ptsrc.size(); ++i)
        result.append("\n"+m_ptsrc[i]->print());

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
 *
 * @todo Do we need the handoff stuff there???
 ***************************************************************************/
void GLATResponse::init_members(void)
{
    // Initialise members
    m_caldb.clear();
    m_rspname.clear();
    m_hasfront   = false;
    m_hasback    = false;
    m_force_mean = false;
    m_aeff.clear();
    m_psf.clear();
    m_edisp.clear();
    m_ptsrc.clear();
    
    // By default use HANDOFF response database.
    char* handoff = std::getenv("HANDOFF_IRF_DIR");
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
    m_caldb      = rsp.m_caldb;
    m_rspname    = rsp.m_rspname;
    m_hasfront   = rsp.m_hasfront;
    m_hasback    = rsp.m_hasback;
    m_force_mean = rsp.m_force_mean;
    m_aeff       = rsp.m_aeff;
    m_psf        = rsp.m_psf;
    m_edisp      = rsp.m_edisp;
    m_ptsrc      = rsp.m_ptsrc;

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


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/
