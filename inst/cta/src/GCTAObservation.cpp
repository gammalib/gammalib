/***************************************************************************
 *               GCTAObservation.cpp  -  CTA Observation class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
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
 * @file GCTAObservation.cpp
 * @brief CTA observation class implementation
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GObservationRegistry.hpp"
#include "GCTAException.hpp"
#include "GCTAObservation.hpp"
#include "GCTAEventList.hpp"
#include "GCTAEventCube.hpp"
#include "GCTARoi.hpp"
#include "GException.hpp"
#include "GFits.hpp"
#include "GTools.hpp"
#include "GModelSpatialPtsrc.hpp"
#include "GIntegral.hpp"
#include "GIntegrand.hpp"

/* __ Globals ____________________________________________________________ */
const GCTAObservation      g_obs_cta_seed;
const GObservationRegistry g_obs_cta_registry(&g_obs_cta_seed);

/* __ Method name definitions ____________________________________________ */
#define G_READ_DS_EBOUNDS       "GCTAObservation::read_ds_ebounds(GFitsHDU*)"
#define G_READ_DS_ROI               "GCTAObservation::read_ds_roi(GFitsHDU*)"

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
 * Creates empty class instance.
 ***************************************************************************/
GCTAObservation::GCTAObservation(void) : GObservation()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] obs CTA observation.
 *
 * Creates class instance by copying an existing CTA observation.
 ***************************************************************************/
GCTAObservation::GCTAObservation(const GCTAObservation& obs) : GObservation(obs)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(obs);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCTAObservation::~GCTAObservation(void)
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
 * @param[in] obs CTA observation.
 *
 * Assign CTA observation to this object.
 ***************************************************************************/
GCTAObservation& GCTAObservation::operator= (const GCTAObservation& obs)
{
    // Execute only if object is not identical
    if (this != &obs) {

        // Copy base class members
        this->GObservation::operator=(obs);

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(obs);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                              Public methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear instance
 ***************************************************************************/
void GCTAObservation::clear(void)
{
    // Free members
    free_members();
    this->GObservation::free_members();

    // Initialise members
    this->GObservation::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
***************************************************************************/
GCTAObservation* GCTAObservation::clone(void) const
{
    return new GCTAObservation(*this);
}


/***********************************************************************//**
 * @brief Set CTA response function
 *
 * @param[in] irfname Name of CTA response function.
 * @param[in] caldb Optional name of calibration database.
 ***************************************************************************/
void GCTAObservation::response(const std::string& irfname, std::string caldb)
{
    // Delete old response function
    if (m_response != NULL) delete m_response;

    // Allocate new CTA response function
    m_response = new GCTAResponse;

    // Set calibration database
    m_response->caldb(caldb);

    // Load instrument response function
    m_response->load(irfname);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns pointer to CTA response function
 ***************************************************************************/
GCTAResponse* GCTAObservation::response(void) const
{
    // Return response pointer
    return m_response;
}


/***********************************************************************//**
 * @brief Returns pointer to CTA pointing direction
 *
 * @param[in] time Time.
 *
 * Returns pointer to pointing direction for a given time.
 *
 * @todo Update pointing information as function of time.
 ***************************************************************************/
GCTAPointing* GCTAObservation::pointing(const GTime& time) const
{
    // Return pointing pointer
    return m_pointing;
}


/***********************************************************************//**
 * @brief Set CTA pointing direction
 *
 * @param[in] pointing Pointing.
 ***************************************************************************/
void GCTAObservation::pointing(const GCTAPointing& pointing)
{
    // Free any existing pointing
    if (m_pointing != NULL) delete m_pointing;
    m_pointing = NULL;

    // Clone pointing
    m_pointing = pointing.clone();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns instrument name
 ***************************************************************************/
std::string GCTAObservation::instrument(void) const
{
    // Return instument name
    return ("CTA");
}


/***********************************************************************//**
 * @brief Print CTA observation information
 ***************************************************************************/
std::string GCTAObservation::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GCTAObservation ===");
    result.append("\n"+parformat("Name")+name());
    result.append("\n"+parformat("Instrument")+instrument());
    result.append("\n"+parformat("Statistics")+statistics());

    // Append pointing
    if (m_pointing != NULL)
        result.append("\n"+m_pointing->print());
    else
        result.append("\n"+parformat("CTA pointing")+"undefined");

    // Append response
    if (m_response != NULL)
        result.append("\n"+response()->print());
    else
        result.append("\n"+parformat("CTA response")+"undefined");

    // Append events
    if (m_events != NULL)
        result.append("\n"+m_events->print());

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Load data for unbinned analysis
 *
 * @param[in] filename Event FITS file name.
 ***************************************************************************/
void GCTAObservation::load_unbinned(const std::string& filename)
{
    // Delete any existing event container (do not call clear() as we do not
    // want to delete the response function)
    if (m_events != NULL) delete m_events;
    m_events = NULL;

    // Allocate event list
    GCTAEventList* events = new GCTAEventList;

    // Assign event list as the observation's event container
    m_events = events;

    // Open FITS file
    GFits file(filename);

    // Read event list
    events->read(file);

    // Read observation attributes from EVENTS extension
    GFitsHDU* hdu = file.hdu("EVENTS");
    read_attributes(hdu);

    // Close FITS file
    file.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load data for binned analysis
 *
 * @param[in] filename Counts map FITS file name.
 ***************************************************************************/
void GCTAObservation::load_binned(const std::string& filename)
{
    // Delete any existing event container (do not call clear() as we do not
    // want to delete the response function)
    if (m_events != NULL) delete m_events;
    m_events = NULL;

    // Allocate event cube
    GCTAEventCube* events = new GCTAEventCube;

    // Assign event cube as the observation's event container
    m_events = events;

    // Open FITS file
    GFits file(filename);

    // Read event cube
    events->read(file);

    // Read observation attributes from primary extension
    GFitsHDU* hdu = file.hdu(0);
    read_attributes(hdu);

    // Close FITS file
    file.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save CTA observation into FITS file.
 *
 * @param[in] filename FITS filename.
 * @param[in] clobber Overwrite existing FITS file (default=false).
 ***************************************************************************/
void GCTAObservation::save(const std::string& filename, bool clobber) const
{
    // Create FITS file
    GFits fits;

    // Get pointers on event list
    GCTAEventList* list = dynamic_cast<GCTAEventList*>(m_events);
    GCTAEventCube* cube = dynamic_cast<GCTAEventCube*>(m_events);

    // Case A: Observation contains an event list
    if (list != NULL) {

        // Write event list into FITS file. This method also writes
        // the GTI as they are part of the event list.
        list->write(fits);

        // Write observation attributes into EVENTS header
        GFitsHDU* hdu = fits.hdu("EVENTS");
        write_attributes(hdu);

    } // endif: observation contained an event list

    // Case B: Observation contains an event cube
    else if (cube != NULL) {

        // Write events cube into FITS file. This method also writes
        // the energy boundaries and the GTI as they are also part
        // of the event cube.
        cube->write(fits);

        // Write observation attributes into primary header
        GFitsHDU* hdu = fits.hdu(0);
        write_attributes(hdu);

    } // endelse: observation contained an event cube

    // Save FITS file
    fits.saveto(filename, clobber);

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
void GCTAObservation::init_members(void)
{
    // Initialise members
    m_response = NULL;
    m_pointing = NULL;
    m_obs_id   = 0;
    m_livetime = 0.0;
    m_ra_obj   = 0.0;
    m_dec_obj  = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] obs CTA observation.
 ***************************************************************************/
void GCTAObservation::copy_members(const GCTAObservation& obs)
{
    // Clone members
    if (obs.m_response != NULL) m_response = obs.m_response->clone();
    if (obs.m_pointing != NULL) m_pointing = obs.m_pointing->clone();

    // Copy members
    m_obs_id   = obs.m_obs_id;
    m_livetime = obs.m_livetime;
    m_ra_obj   = obs.m_ra_obj;
    m_dec_obj  = obs.m_dec_obj;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAObservation::free_members(void)
{
    // Free memory
    if (m_response != NULL) delete m_response;
    if (m_pointing != NULL) delete m_pointing;

    // Mark memory as free
    m_response = NULL;
    m_pointing = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read observation attributes
 *
 * @param[in] hdu FITS HDU pointer
 *
 * Reads the observation attributes from HDU. Nothing is done if the HDU
 * pointer is NULL.
 ***************************************************************************/
void GCTAObservation::read_attributes(const GFitsHDU* hdu)
{
    // Continue only if HDU is valid
    if (hdu != NULL) {

        // Read attributes
        m_obs_id   = hdu->integer("OBS_ID");
        m_livetime = hdu->real("LIVETIME");
        m_name     = hdu->string("OBJECT");
        m_ra_obj   = hdu->real("RA_OBJ");
        m_dec_obj  = hdu->real("DEC_OBJ");

        // Read pointing information
        double  ra_pnt  = hdu->real("RA_PNT");
        double  dec_pnt = hdu->real("DEC_PNT");
        double  alt_pnt = hdu->real("ALT_PNT");
        double  az_pnt  = hdu->real("AZ_PNT");
        GSkyDir pnt;
        pnt.radec_deg(ra_pnt, dec_pnt);

        // Set pointing
        if (m_pointing != NULL) delete m_pointing;
        m_pointing = new GCTAPointing;
        m_pointing->dir(pnt);

    } // endif: HDU was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write observation attributes
 *
 * @param[in] hdu FITS HDU pointer
 *
 * Nothing is done if the HDU pointer is NULL.
 ***************************************************************************/
void GCTAObservation::write_attributes(GFitsHDU* hdu) const
{
    // Continue only if HDU is valid
    if (hdu != NULL) {

        // Compute some attributes
        double ra_pnt  = (m_pointing != NULL) ? m_pointing->dir().ra_deg() : 0.0;
        double dec_pnt = (m_pointing != NULL) ? m_pointing->dir().dec_deg() : 0.0;
        double tstart  = events()->tstart().met();
        double tstop   = events()->tstop().met();
        double telapse = events()->gti().telapse();
        double ontime  = events()->gti().ontime();
        double deadc   = (ontime > 0.0) ? livetime() / ontime : 0.0;

        // Set observation information
        hdu->card("CREATOR",  "GammaLib", "Program which created the file");
        hdu->card("TELESCOP", "CTA",      "Telescope");
        hdu->card("OBS_ID",   obs_id(),   "Observation identifier");
        hdu->card("DATE_OBS", "string",   "Observation start date");
        hdu->card("TIME_OBS", "string",   "Observation start time");
        hdu->card("DATE_END", "string",   "Observation end date");
        hdu->card("TIME_END", "string",   "Observation end time");

        // Set observation time information
        hdu->card("TSTART",   tstart, "[s] Mission time of start of observation");
        hdu->card("TSTOP",    tstop, "[s] Mission time of end of observation");
        hdu->card("MJDREFI",  51910, "[days] Integer part of mission time reference MJD");
        hdu->card("MJDREFF",  7.428703703703703e-14, "[days] Fractional part of mission time reference MJD");
        hdu->card("TIMEUNIT", "s", "Time unit");
        hdu->card("TIMESYS",  "TT", "Time system");
        hdu->card("TIMEREF",  "LOCAL", "Time reference");
        hdu->card("TELAPSE",  telapse, "[s] Mission elapsed time");
        hdu->card("ONTIME",   ontime, "[s] Total good time including deadtime");
        hdu->card("LIVETIME", livetime(), "[s] Total livetime");
        hdu->card("DEADC",    deadc, "Deadtime fraction");
        hdu->card("TIMEDEL",  1.0, "Time resolution");

        // Set pointing information
        hdu->card("OBJECT",   name(),    "Observed object");
        hdu->card("RA_OBJ",   ra_obj(),  "[deg] Target Right Ascension");
        hdu->card("DEC_OBJ",  dec_obj(), "[deg] Target Declination");
        hdu->card("RA_PNT",   ra_pnt,    "[deg] Pointing Right Ascension");
        hdu->card("DEC_PNT",  dec_pnt,   "[deg] Pointing Declination");
        hdu->card("ALT_PNT",  0.0,       "[deg] Average altitude of pointing");
        hdu->card("AZ_PNT",   0.0,       "[deg] Average azimuth of pointing");
        hdu->card("RADECSYS", "FK5",     "Coordinate system");
        hdu->card("EQUINOX",  2000.0,    "Epoch");
        hdu->card("CONV_DEP", 0.0,       "Convergence depth of telescopes");
        hdu->card("CONV_RA",  0.0,       "[deg] Convergence Right Ascension");
        hdu->card("CONV_DEC", 0.0,       "[deg] Convergence Declination");
        hdu->card("OBSERVER", "string",  "Observer");

        // Telescope information
        hdu->card("N_TELS",   100,      "Number of telescopes in event list");
        hdu->card("TELLIST",  "string", "Telescope IDs");
        hdu->card("GEOLAT",   0.0,      "[deg] Geographic latitude of array centre");
        hdu->card("GEOLON",   0.0,      "[deg] Geographic longitude of array centre");
        hdu->card("ALTITUDE", 0.0,      "[km] Altitude of array centre");

        // Other information
        hdu->card("EUNIT",    "TeV",    "Energy unit");
        hdu->card("EVTVER",   "draft1", "Event list version number");

    } // endif: HDU was valid

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                        Npred integration methods                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Temporally integrate spatially & spectrally integrated Npred kernel
 *
 * @param[in] model Gamma-ray source model.
 *
 * Implement the temporal integration as a simple multiplication by the
 * elapsed time. This assumes that the source is non-variable during the
 * observation and that the CTA pointing is stable.
 ***************************************************************************/
double GCTAObservation::npred_temp(const GModel& model) const
{
    // Initialise result
    double result = 0.0;

    // Determine ontime
    double ontime = events()->gti().ontime();

    // Integrate only if ontime is positive
    if (ontime > 0.0) {

        // Integration is a simple multiplication by the time
        result = npred_spec(model, events()->gti().tstart()) * ontime;

    }

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/
