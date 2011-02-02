/***************************************************************************
 *               GCTAObservation.cpp  -  CTA Observation class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GCTAObservation.cpp
 * @brief GCTAObservation class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
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
 * @brief Constructor
 *
 * Creates an empty instance of GCTAObservation.
 ***************************************************************************/
GCTAObservation::GCTAObservation(void) : GObservation()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] obs Observation from which the instance should be built.
 *
 * Creates an instance of GCTAObservation by copying information from an
 * existing instance.
 ***************************************************************************/
GCTAObservation::GCTAObservation(const GCTAObservation& obs) : GObservation(obs)
{
    // Initialise class members for clean destruction
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
 * @param[in] obs Observation which should be assigned.
 *
 * Assign instance of GCTAObservation to this object.
 ***************************************************************************/
GCTAObservation& GCTAObservation::operator= (const GCTAObservation& obs)
{
    // Execute only if object is not identical
    if (this != &obs) {

        // Copy base class members
        this->GObservation::operator=(obs);

        // Free members
        free_members();

        // Initialise private members for clean destruction
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
    result.append("=== GCTAObservation ===\n");
    result.append(parformat("Name")+obsname()+"\n");
    result.append(parformat("Instrument")+instrument()+"\n");
    result.append(parformat("Statistics")+statistics()+"\n");

    // Append time range
    result.append(parformat("Time range"));
    result.append(str(m_gti.tstart().met()));
    result.append(" - ");
    result.append(str(m_gti.tstop().met()));
    result.append(" s\n");

    // Append energy range
    result.append(parformat("Energy range"));
    result.append(m_ebounds.emin().print());
    result.append(" - ");
    result.append(m_ebounds.emax().print());

    // Append pointing
    if (m_pointing != NULL)
        result.append("\n"+m_pointing->print());
    else
        result.append("\n"+parformat("CTA pointing")+"undefined");

    // Append ROI
    if (m_roi != NULL)
        result.append("\n"+roi()->print());
    else
        result.append("\n"+parformat("Region of interest")+"undefined");

    // Append GTIs
    //result.append("\n"+gti().print());

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
    // Delete old events and ROI. We do not call clear() here since we want
    // to preserve any existing response function.
    if (m_events != NULL) delete m_events;
    if (m_roi    != NULL) delete m_roi;
    m_events = NULL;
    m_roi    = NULL;

    // Reset energy boundaries
    m_ebounds.clear();

    // Allocate events
    GCTAEventList* events = new GCTAEventList;
    m_events = events;

    // Open FITS file
    GFits file(filename);

    // Get HDUs
    GFitsTable* hdu = file.table("EVENTS");

    // Read events into list
    events->read(hdu);

    // Read observation attributes
    read_attributes(hdu);

    // Read data selection keywords
    read_ds_ebounds(hdu);
    read_ds_roi(hdu);

    // Close FITS file
    file.close();

    // Load GTIs
    m_gti.load(filename);

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
    // Delete old events. We do not call clear() here since we want to
    // preserve any existing response function.
    if (m_events != NULL) delete m_events;

    // Allocate events
    GCTAEventCube* events = new GCTAEventCube;
    m_events = events;

    // Read events into list
    events->load(filename);

    // Read observation attributes from first header
    GFits file(filename);
    GFitsHDU* hdu = file.hdu(0);
    read_attributes(hdu);
    file.close();

    // Copy energy boundaries and GTIs from event cube
    m_ebounds = events->ebounds();
    m_gti     = events->gti();

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

    }

    // Case B: Observation contains an event cube
    else if (cube != NULL) {

        // Copy energy boundaries and GTIs into event cube
        cube->ebounds(m_ebounds);
        cube->gti(m_gti);

        // Write events cube into FITS file
        cube->write(&fits);

        // Write observation attributes into first header
        GFitsHDU* hdu = fits.hdu(0);
        write_attributes(hdu);

    } // endif: observation contained an events cube

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
 *
 * @todo We allocate void response and pointing instances so make sure they
 * exist (analysis methods depend on the existence of these members).
 ***************************************************************************/
void GCTAObservation::init_members(void)
{
    // Initialise members
    m_response = NULL;
    m_pointing = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] obs Observation to be copied
 *
 * @todo Try to avoid the back pointer if possible, or flesh out a way that
 * makes this more solid. The point is: when events are copied the back
 * pointer to the observation needs to be updated, yet when the copying is
 * done in the events class the class does not known about the observation.
 * Thus, back pointer update has to be done by the observation class.
 ***************************************************************************/
void GCTAObservation::copy_members(const GCTAObservation& obs)
{
    // Copy members
    if (obs.m_response != NULL) m_response = obs.m_response->clone();
    if (obs.m_pointing != NULL) m_pointing = obs.m_pointing->clone();

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
 * @param[in] hdu FITS HDU
 *
 * Reads the observation attributes.
 ***************************************************************************/
void GCTAObservation::read_attributes(GFitsHDU* hdu)
{
    // Continue only if HDU is valid
    if (hdu != NULL) {

        // Read attributes
        m_obs_id   = hdu->integer("OBS_ID");
        m_livetime = hdu->real("LIVETIME");
        m_obsname  = hdu->string("OBJECT");
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
 * @brief Read energy boundary data selection keywords
 *
 * @param[in] hdu FITS HDU
 *
 * @exception GCTAException::no_ebds
 *            Invalid energy data selection encountered
 *
 * Reads the energy boundary data selection keywords by searching for a 
 * DSTYPx keyword named "ENERGY". The data selection information is expected
 * to be in the format "200:50000", where the 2 arguments are the minimum
 * and maximum energy. The energy unit is given by the keyword DSUNIx, which
 * supports keV, MeV, GeV and TeV (case independent). No detailed syntax
 * checking is performed.
 *
 * This method does nothing if the HDU is NULL.
 ***************************************************************************/
void GCTAObservation::read_ds_ebounds(GFitsHDU* hdu)
{
    // Continue only if HDU is valid
    if (hdu != NULL) {

        // Reset energy boundaries
        m_ebounds.clear();

        // Get number of data selection keywords
        int ndskeys = 0;
        try {
            ndskeys = hdu->integer("NDSKEYS");
        }
        catch (GException::fits_key_not_found) {
            ;
        }

        // Loop over all data selection keys
        for (int i = 1; i <= ndskeys; ++i) {
            std::string type_key  = "DSTYP"+str(i);
            std::string unit_key  = "DSUNI"+str(i);
            std::string value_key = "DSVAL"+str(i);
            try {
                if (hdu->string(type_key) == "ENERGY") {
                    std::string unit                 = toupper(hdu->string(unit_key));
                    std::string value                = hdu->string(value_key);
                    std::vector<std::string> ebounds = split(value, ":");
                    if (ebounds.size() == 2) {
                        double  emin = todouble(ebounds[0]);
                        double  emax = todouble(ebounds[1]);
                        GEnergy e_min;
                        GEnergy e_max;
                        if (unit == "KEV") {
                            e_min.keV(emin);
                            e_max.keV(emax);
                        }
                        else if (unit == "MEV") {
                            e_min.MeV(emin);
                            e_max.MeV(emax);
                        }
                        else if (unit == "GEV") {
                            e_min.GeV(emin);
                            e_max.GeV(emax);
                        }
                        else if (unit == "TEV") {
                            e_min.TeV(emin);
                            e_max.TeV(emax);
                        }
                        else {
                            throw GCTAException::no_ebds(G_READ_DS_EBOUNDS,
                                  "Invalid energy unit \""+unit+
                                  "\" encountered in data selection key \""+
                                  unit_key+"\"");
                        }
                        m_ebounds.append(e_min, e_max);
                    }
                    else {
                        throw GCTAException::no_ebds(G_READ_DS_EBOUNDS,
                              "Invalid energy value \""+value+
                              "\" encountered in data selection key \""+
                              value_key+"\"");
                    }
                } // endif: ENERGY type found
            }
            catch (GException::fits_key_not_found) {
                ;
            }
        } // endfor: looped over data selection keys

    } // endif: HDU was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read ROI data selection keywords
 *
 * @param[in] hdu FITS HDU
 *
 * @exception GException::no_roi
 *            Invalid ROI data selection encountered
 *
 * Reads the ROI data selection keywords by searching for a DSTYPx keyword
 * named "POS(RA,DEC)". The data selection information is expected to be
 * in the format "CIRCLE(267.0208,-24.78,4.5)", where the 3 arguments are
 * Right Ascension, Declination and radius in units of degrees. No detailed
 * syntax checking is performed.
 *
 * This method does nothing if the HDU is NULL.
 ***************************************************************************/
void GCTAObservation::read_ds_roi(GFitsHDU* hdu)
{
    // Continue only if HDU is valid
    if (hdu != NULL) {

        // Delete any old ROI
        if (m_roi    != NULL) delete m_roi;

        // Allocate ROI
        GCTARoi* roi = new GCTARoi;
        m_roi = roi;

        // Get number of data selection keywords
        int ndskeys = 0;
        try {
            ndskeys = hdu->integer("NDSKEYS");
        }
        catch (GException::fits_key_not_found) {
            ;
        }

        // Loop over all data selection keys
        for (int i = 1; i <= ndskeys; ++i) {
            std::string type_key  = "DSTYP"+str(i);
            std::string unit_key  = "DSUNI"+str(i);
            std::string value_key = "DSVAL"+str(i);
            try {
                if (hdu->string(type_key) == "POS(RA,DEC)") {
                    std::string unit              = toupper(hdu->string(unit_key));
                    std::string value             = hdu->string(value_key);
                    value                         = strip_chars(value, "CIRCLE(");
                    value                         = strip_chars(value, ")");
                    std::vector<std::string> args = split(value, ",");
                    if (args.size() == 3) {
                        double  ra  = todouble(args[0]);
                        double  dec = todouble(args[1]);
                        double  rad = todouble(args[2]);
                        GCTAInstDir dir;
                        dir.radec_deg(ra, dec);
                        roi->centre(dir);
                        roi->radius(rad);
                    }
                    else {
                        throw GException::no_roi(G_READ_DS_ROI,
                              "Invalid acceptance cone value \""+value+
                              "\" encountered in data selection key \""+
                              value_key+"\"");
                    }
                } // endif: POS(RA,DEC) type found
            }
            catch (GException::fits_key_not_found) {
                ;
            }
        } // endfor: looped over data selection keys

    } // endif: HDU was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write observation attributes
 *
 * @param[in] hdu FITS HDU
 ***************************************************************************/
void GCTAObservation::write_attributes(GFitsHDU* hdu) const
{
    // Continue only if HDU is valid
    if (hdu != NULL) {

        // Compute some attributes
        double deadc   = livetime() / gti().ontime();
        double ra_pnt  = (m_pointing != NULL) ? m_pointing->dir().ra_deg() : 0.0;
        double dec_pnt = (m_pointing != NULL) ? m_pointing->dir().dec_deg() : 0.0;
        double tstart  = this->tstart().met();
        double tstop   = this->tstop().met();

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
        hdu->card("TELAPSE",  gti().telapse(), "[s] Mission elapsed time");
        hdu->card("ONTIME",   gti().ontime(), "[s] Total good time including deadtime");
        hdu->card("LIVETIME", livetime(), "[s] Total livetime");
        hdu->card("DEADC",    deadc, "Deadtime fraction");
        hdu->card("TIMEDEL",  1.0, "Time resolution");

        // Set pointing information
        hdu->card("OBJECT",   obsname(), "Observed object");
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
    double ontime = m_gti.ontime();

    // Integrate only if ontime is positive
    if (ontime > 0.0) {

        // Integration is a simple multiplication by the time
        result = npred_spec(model, m_gti.tstart()) * ontime;

    }

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/
