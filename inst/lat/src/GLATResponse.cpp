/***************************************************************************
 *                GLATResponse.cpp - Fermi LAT response class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2016 by Juergen Knoedlseder                         *
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
 * @file GLATResponse.cpp
 * @brief Fermi LAT response class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <unistd.h>           // access() function
#include <cstdlib>            // std::getenv() function
#include <string>
#include "GException.hpp"
#include "GFits.hpp"
#include "GTools.hpp"
#include "GCaldb.hpp"
#include "GSource.hpp"
#include "GModelSpatialPointSource.hpp"
#include "GLATInstDir.hpp"
#include "GLATResponse.hpp"
#include "GLATObservation.hpp"
#include "GLATEventAtom.hpp"
#include "GLATEventBin.hpp"
#include "GLATEventCube.hpp"
#include "GLATException.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_CALDB                           "GLATResponse::caldb(std::string&)"
#define G_LOAD                             "GLATResponse::load(std::string&)"
#define G_IRF      "GLATResponse::irf(GInstDir&, GEnergy&, GTime&, GSkyDir&,"\
                                          " GEnergy&, GTime&, GObservation&)"
#define G_NROI            "GLATResponse::nroi(GModelSky&, GEnergy&, GTime&, "\
                                                             "GObservation&)"
#define G_EBOUNDS                           "GLATResponse::ebounds(GEnergy&)"
#define G_AEFF                                     "GLATResponse::aeff(int&)"
#define G_PSF                                       "GLATResponse::psf(int&)"
#define G_EDISP                                   "GLATResponse::edisp(int&)"
#define G_IRF_ATOM     "GLATResponse::irf(GLATEventAtom&, GModel&, GEnergy&,"\
                                                     "GTime&, GObservation&)"
#define G_IRF_BIN       "GLATResponse::irf(GLATEventBin&, GModel&, GEnergy&,"\
                                                     "GTime&, GObservation&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
#define G_DUMP_MEAN_PSF  0                    //!< Dump mean PSF allocation
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
    // Initialise members
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
    // Initialise members
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
 * @return Response.
 ***************************************************************************/
GLATResponse& GLATResponse::operator=(const GLATResponse& rsp)
{
    // Execute only if object is not identical
    if (this != &rsp) {

        // Copy base class members
        this->GResponse::operator=(rsp);

        // Free members
        free_members();

        // Initialise members
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
 * @brief Clear response
***************************************************************************/
void GLATResponse::clear(void)
{
    // Free class members
    free_members();
    this->GResponse::free_members();

    // Initialise members
    this->GResponse::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone response
 *
 * @return Pointer to deep copy of response.
 ***************************************************************************/
GLATResponse* GLATResponse::clone(void) const
{
    return new GLATResponse(*this);
}


/***********************************************************************//**
 * @brief Return value of point source IRF
 *
 * @param[in] event Observed event.
 * @param[in] photon Incident photon.
 * @param[in] obs Observation.
 *
 * @exception GLATException::bad_instdir_type
 *            Instrument direction is not a valid LAT instrument direction.
 *
 * @todo The IRF value is not divided by ontime of the event, but it is
 *       already time integrated.
 ***************************************************************************/
double GLATResponse::irf(const GEvent&       event,
                         const GPhoton&      photon,
                         const GObservation& obs) const
{
    // Get pointer to LAT instrument direction
    const GLATInstDir* dir = dynamic_cast<const GLATInstDir*>(&(event.dir()));
    if (dir == NULL) {
        throw GLATException::bad_instdir_type(G_IRF);
    }

    // Get photon attributes
    const GSkyDir& srcDir = photon.dir();
    const GEnergy& srcEng = photon.energy();

    // Search for mean PSF
    int ipsf = -1;
    for (int i = 0; i < m_ptsrc.size(); ++i) {
        if (m_ptsrc[i]->dir() == srcDir) {
            ipsf = i;
            break;
        }
    }

    // If mean PSF has not been found then create it now
    if (ipsf == -1) {

        // Allocate new mean PSF
        GLATMeanPsf* psf = new GLATMeanPsf(srcDir, static_cast<const GLATObservation&>(obs));

        // Set source name
        std::string name = "SRC("+gammalib::str(srcDir.ra_deg()) +
                                "," +
                                gammalib::str(srcDir.dec_deg())+")";
        psf->name(name);

        // Push mean PSF on stack
        const_cast<GLATResponse*>(this)->m_ptsrc.push_back(psf);

        // Set index of mean PSF
        ipsf = m_ptsrc.size()-1;

        // Debug option: dump mean PSF
        #if G_DUMP_MEAN_PSF
        std::cout << "Added new mean PSF \""+name+"\"" << std::endl;
        std::cout << *psf << std::endl;
        #endif

    } // endif: created new mean PSF

    // Get IRF value
    double offset = dir->dir().dist_deg(srcDir);
    double irf    = (*m_ptsrc[ipsf])(offset, srcEng.log10MeV());

    // Return IRF value
    return irf;
}


/***********************************************************************//**
 * @brief Return value of model instrument response function
 *
 * @param[in] event Event.
 * @param[in] source Source.
 * @param[in] obs Observation.
 *
 * This method returns the response of the instrument to a specific source
 * model. The method handles both event atoms and event bins.
 ***************************************************************************/
double GLATResponse::irf(const GEvent&       event,
                         const GSource&      source,
                         const GObservation& obs) const
{
    // Get IRF value
    double rsp;
    if (event.is_atom()) {
        rsp = irf(static_cast<const GLATEventAtom&>(event), source, obs);
    }
    else {
        rsp = irf(static_cast<const GLATEventBin&>(event), source, obs);
    }

    // Return IRF value
    return rsp;
}


/***********************************************************************//**
 * @brief Return value of model IRF for event atom
 *
 * @param[in] event Event atom.
 * @param[in] source Source.
 * @param[in] obs Observations.
 *
 * @exception GException::feature_not_implemented
 *            Method not yet implemented.
 *
 * @todo Not yet implemented.
 ***************************************************************************/
double GLATResponse::irf(const GLATEventAtom& event,
                         const GSource&       source,
                         const GObservation&  obs) const
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
 * @param[in] source Source.
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
double GLATResponse::irf(const GLATEventBin& event,
                         const GSource&      source,
                         const GObservation& obs) const
{
    // Initialise response value
    double rsp = 0.0;

    // Get pointer to event cube
    const GLATEventCube* cube = event.cube();

    // Get pointer on point source spatial model
    const GModelSpatialPointSource* ptsrc =
          dynamic_cast<const GModelSpatialPointSource*>(source.model());

    // Get source energy
    GEnergy srcEng = source.energy();

    // Search for diffuse response in event cube
    int idiff = -1;
    for (int i = 0; i < cube->ndiffrsp(); ++i) {
        if (cube->diffname(i) == source.name()) {
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
        GSkyMap*      map    = cube->diffrsp(idiff);
        const double* pixels = map->pixels() + event.ipix();
        rsp                  = nodes.wgt_left()  * pixels[nodes.inx_left()  * map->npix()] +
                               nodes.wgt_right() * pixels[nodes.inx_right() * map->npix()];

        // Divide by solid angle and ontime since source maps are given in units of
        // counts/pixel/MeV.
        rsp /= (event.solidangle() * event.ontime());

    } // endelse: diffuse response was present

    // ... otherwise check if model is a point source. If this is true
    // then return response from mean PSF
    if ((idiff == -1 || m_force_mean) && ptsrc != NULL) {

        // Search for mean PSF
        int ipsf = -1;
        for (int i = 0; i < m_ptsrc.size(); ++i) {
            if (m_ptsrc[i]->name() == source.name()) {
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
            psf->name(source.name());

            // Push mean PSF on stack
            const_cast<GLATResponse*>(this)->m_ptsrc.push_back(psf);

            // Set index of mean PSF
            ipsf = m_ptsrc.size()-1;

            // Debug option: dump mean PSF
            #if G_DUMP_MEAN_PSF 
            std::cout << "Added new mean PSF \""+source.name() << "\"" << std::endl;
            std::cout << *psf << std::endl;
            #endif

        } // endif: created new mean PSF

        // Get PSF value
        GSkyDir srcDir   = m_ptsrc[ipsf]->dir();
        double  offset   = event.dir().dir().dist_deg(srcDir);
        double  mean_psf = (*m_ptsrc[ipsf])(offset, srcEng.log10MeV()) / (event.ontime());

        // Debug option: compare mean PSF to diffuse response
        #if G_DEBUG_MEAN_PSF
        std::cout << "Energy=" << srcEng.MeV();
        std::cout << " MeanPsf=" << mean_psf;
        std::cout << " DiffusePsf=" << rsp;
        std::cout << " ratio(Mean/Diffuse)=" << mean_psf/rsp;
        if (mean_psf/rsp < 0.99) {
            std::cout << " <<< (1%)";
        }
        else if (mean_psf/rsp > 1.01) {
            std::cout << " >>> (1%)";
        }
        std::cout << std::endl;
        #endif

        // Set response
        rsp = mean_psf;

    } // endif: model was point source

    // ... otherwise throw an exception
    if ((idiff == -1) && ptsrc == NULL) {
        throw GLATException::diffuse_not_found(G_IRF_BIN, source.name());
    }

    // Return IRF value
    return rsp;
}


/***********************************************************************//**
 * @brief Return integral of event probability for a given sky model over ROI
 *
 * @param[in] model Sky model.
 * @param[in] obsEng Observed photon energy.
 * @param[in] obsTime Observed photon arrival time.
 * @param[in] obs Observation.
 * @return 0
 *
 * @exception GException::feature_not_implemented
 *            Method is not implemented.
 ***************************************************************************/
double GLATResponse::nroi(const GModelSky&    model,
                          const GEnergy&      obsEng,
                          const GTime&        obsTime,
                          const GObservation& obs) const
{
    // Method is not implemented
    std::string msg = "Spatial integration of sky model over the data space "
                      "is not implemented.";
    throw GException::feature_not_implemented(G_NROI, msg);

    // Return Npred
    return (0.0);
}


/***********************************************************************//**
 * @brief Return true energy boundaries for a specific observed energy
 *
 * @param[in] obsEnergy Observed Energy.
 * @return True energy boundaries for given observed energy.
 *
 * @exception GException::feature_not_implemented
 *            Method is not implemented.
 ***************************************************************************/
GEbounds GLATResponse::ebounds(const GEnergy& obsEnergy) const
{
    // Initialise an empty boundary object
    GEbounds ebounds;

    // Throw an exception
    std::string msg = "Energy dispersion not implemented.";
    throw GException::feature_not_implemented(G_EBOUNDS, msg);

    // Return energy boundaries
    return (ebounds);
}


/***********************************************************************//**
 * @brief Load Fermi LAT response from calibration database
 *
 * @param[in] rspname Response name.
 *
 * @exception GException::invalid_argument
 *            Invalid response type encountered.
 *
 * Loads the specified Fermi LAT response from the calibration database.
 * The following response names are supported (case insensitive):
 *
 *      name (is equivalent to front+back)
 *      name::front
 *      name::back
 *      name::psf(0-3)
 *      name::edisp(0-3)
 *
 * where name is the response name (for example "P8R2_SOURCE_V6"). Note that
 * the name is case sensitive, but the event typ is not case sensitive.
 ***************************************************************************/
void GLATResponse::load(const std::string& rspname)
{
    // Set possible response types
    static const std::string type_names[] = {"FRONT",
                                             "BACK",
                                             "PSF0",
                                             "PSF1",
                                             "PSF2",
                                             "PSF3",
                                             "EDISP0",
                                             "EDISP1",
                                             "EDISP2",
                                             "EDISP3"};

    // Clear instance but conserve calibration database
    GCaldb caldb = m_caldb;
    clear();
    m_caldb = caldb;

    // Determine response types to be loaded. If no event type has been
    // specified then load both front and back response.
    int event_type = 0;
    std::vector<std::string> array = gammalib::split(rspname, ":");
    if (array.size() == 1) {
        m_rspname   = array[0];
        event_type  = 3;      // 0x0000000011
    }
    else if (array.size() == 3) {
        m_rspname          = array[0];
        std::string evtype = gammalib::strip_whitespace(gammalib::tolower(array[2]));
        if (evtype == "front") {
            event_type = 1;   // 0x0000000001
        }
        else if (evtype == "back") {
            event_type = 2;   // 0x0000000010
        }
        else if (evtype == "psf0") {
            event_type = 4;   // 0x0000000100
        }
        else if (evtype == "psf1") {
            event_type = 8;   // 0x0000001000
        }
        else if (evtype == "psf2") {
            event_type = 16;  // 0x0000010000
        }
        else if (evtype == "psf3") {
            event_type = 32;  // 0x0000100000
        }
        else if (evtype == "psf") {
            event_type = 60;  // 0x0000111100
        }
        else if (evtype == "edisp0") {
            event_type = 64;  // 0x0001000000
        }
        else if (evtype == "edisp1") {
            event_type = 128; // 0x0010000000
        }
        else if (evtype == "edisp2") {
            event_type = 256; // 0x0100000000
        }
        else if (evtype == "edisp3") {
            event_type = 512; // 0x1000000000
        }
        else if (evtype == "edisp") {
            event_type = 960; // 0x1111000000
        }
        else {
            std::string msg = "Invalid response type \""+array[2]+"\". "
                              "Expect one of \"FRONT\", \"BACK\", "
                              "\"PSF(0-3)\", or \"EDISP(0-3)\" "
                              "(case insensitive).";
            throw GException::invalid_argument(G_LOAD, msg);
        }
    }
    else {
        std::string msg = "Invalid response \""+m_rspname+"\".";
        throw GException::invalid_argument(G_LOAD, msg);
    }

    // Loop over all possible response types and append the relevant
    // response function components to the response
    for (int i = 0; i < 10; ++i) {

        // Set bitmask for this response type
        int bitmask = 1 << i;

        // Fall through if this response type is not requested
        if ((event_type & bitmask) == 0) {
            continue;
        }

        // Get response using the GCaldb interface
        std::string expr      = "VERSION("+gammalib::toupper(m_rspname)+")";
        std::string aeffname  = m_caldb.filename(type_names[i],"","EFF_AREA","","",expr);
        std::string psfname   = m_caldb.filename(type_names[i],"","RPSF","","",expr);
        std::string edispname = m_caldb.filename(type_names[i],"","EDISP","","",expr);

        // Load IRF components
        GLATAeff*  aeff  = new GLATAeff(aeffname, type_names[i]);
        GLATPsf*   psf   = new GLATPsf(psfname, type_names[i]);
        GLATEdisp* edisp = new GLATEdisp(edispname, type_names[i]);

        // Push IRF components on the response stack
        m_aeff.push_back(aeff);
        m_psf.push_back(psf);
        m_edisp.push_back(edisp);

    } // endfor: looped over response types

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save Fermi LAT response in calibration database
 *
 * @param[in] rspname Response name.
 *
 * @todo Not yet implemented.
 ***************************************************************************/
void GLATResponse::save(const std::string& rspname) const
{
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
    if (index < 0 || index >= m_aeff.size()) {
        throw GException::out_of_range(G_AEFF, index, 0, m_aeff.size()-1);
    }
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
    if (index < 0 || index >= m_psf.size()) {
        throw GException::out_of_range(G_PSF, index, 0, m_psf.size()-1);
    }
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
    if (index < 0 || index >= m_edisp.size()) {
        throw GException::out_of_range(G_EDISP, index, 0, m_edisp.size()-1);
    }
    #endif

    // Return pointer on energy dispersion
    return m_edisp[index];
}


/***********************************************************************//**
 * @brief Print Fermi-LAT response information
 *
 * @param[in] chatter Chattiness.
 * @return String containing response information.
 ***************************************************************************/
std::string GLATResponse::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GLATResponse ===");

        // Append information
        result.append("\n"+gammalib::parformat("Caldb mission")+m_caldb.mission());
        result.append("\n"+gammalib::parformat("Caldb instrument")+m_caldb.instrument());
        result.append("\n"+gammalib::parformat("Response name")+m_rspname);
        result.append("\n"+gammalib::parformat("Event types"));
        for (int i = 0; i < size(); ++i) {
            result.append(m_aeff[i]->evtype()+" ");
        }

        if (chatter > TERSE) {
            for (int i = 0; i < size(); ++i) {
                result.append("\n"+m_aeff[i]->print(gammalib::reduce(chatter)));
                result.append("\n"+m_psf[i]->print(gammalib::reduce(chatter)));
                result.append("\n"+m_edisp[i]->print(gammalib::reduce(chatter)));
            }
            for (int i = 0; i < m_ptsrc.size(); ++i) {
                result.append("\n"+m_ptsrc[i]->print(gammalib::reduce(chatter)));
            }
        }

    } // endif: chatter was not silent

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
    m_caldb.clear();
    m_rspname.clear();
    m_force_mean = false;
    m_aeff.clear();
    m_psf.clear();
    m_edisp.clear();
    m_ptsrc.clear();

    // Open Fermi/LAT calibration database
    m_caldb.open("glast", "lat");
    
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
    m_force_mean = rsp.m_force_mean;

    // Clone Aeff
    m_aeff.clear();
    for (int i = 0; i < rsp.m_aeff.size(); ++i) {
        m_aeff.push_back(rsp.m_aeff[i]->clone());
    }
    
    // Clone Psf
    m_psf.clear();
    for (int i = 0; i < rsp.m_psf.size(); ++i) {
        m_psf.push_back(rsp.m_psf[i]->clone());
    }

    // Clone Edisp
    m_edisp.clear();
    for (int i = 0; i < rsp.m_edisp.size(); ++i) {
        m_edisp.push_back(rsp.m_edisp[i]->clone());
    }

    // Clone point sources
    m_ptsrc.clear();
    for (int i = 0; i < rsp.m_ptsrc.size(); ++i) {
        m_ptsrc.push_back(rsp.m_ptsrc[i]->clone());
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GLATResponse::free_members(void)
{
    // Free Aeff memory
    for (int i = 0; i < m_aeff.size(); ++i) {
        if (m_aeff[i] != NULL) delete m_aeff[i];
        m_aeff[i] = NULL;
    }
    m_aeff.clear();

    // Free Psf memory
    for (int i = 0; i < m_psf.size(); ++i) {
        if (m_psf[i] != NULL) delete m_psf[i];
        m_psf[i] = NULL;
    }
    m_psf.clear();

    // Free Edisp memory
    for (int i = 0; i < m_edisp.size(); ++i) {
        if (m_edisp[i] != NULL) delete m_edisp[i];
        m_edisp[i] = NULL;
    }
    m_edisp.clear();

    // Free point source memory
    for (int i = 0; i < m_ptsrc.size(); ++i) {
        if (m_ptsrc[i] != NULL) delete m_ptsrc[i];
        m_ptsrc[i] = NULL;
    }
    m_ptsrc.clear();

    // Return
    return;
}
