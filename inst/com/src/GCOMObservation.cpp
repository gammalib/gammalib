/***************************************************************************
 *             GCOMObservation.cpp - COMPTEL Observation class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2017 by Juergen Knoedlseder                         *
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
 * @file GCOMObservation.cpp
 * @brief COMPTEL observation class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <typeinfo> 
#include "GObservationRegistry.hpp"
#include "GTools.hpp"
#include "GException.hpp"
#include "GFits.hpp"
#include "GCaldb.hpp"
#include "GCOMObservation.hpp"
#include "GCOMEventCube.hpp"
#include "GCOMEventList.hpp"
#include "GCOMStatus.hpp"
#include "GCOMSupport.hpp"

/* __ Globals ____________________________________________________________ */
const GCOMObservation      g_obs_com_seed;
const GObservationRegistry g_obs_com_registry(&g_obs_com_seed);

/* __ Method name definitions ____________________________________________ */
#define G_RESPONSE                    "GCOMObservation::response(GResponse&)"
#define G_READ                          "GCOMObservation::read(GXmlElement&)"
#define G_WRITE                        "GCOMObservation::write(GXmlElement&)"
#define G_COMPUTE_DRE               "GCOMObservation::compute_dre(GEbounds&)"
#define G_LOAD_DRB                    "GCOMObservation::load_drb(GFilename&)"
#define G_LOAD_DRG                    "GCOMObservation::load_drg(GFilename&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
#define G_CHECK_EHA_COMPUTATION

/* __ Debug definitions __________________________________________________ */
#define G_DEBUG_COMPUTE_DRE


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 *
 * Creates an empty COMPTEL observation.
 ***************************************************************************/
GCOMObservation::GCOMObservation(void) : GObservation()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Binned observation constructor
 *
 * @param[in] drename Event cube name.
 * @param[in] drbname Background cube name.
 * @param[in] drgname Geometry cube name.
 * @param[in] drxname Exposure map name.
 *
 * Creates COMPTEL observation by loading the following FITS files:
 *
 *     DRE - Events cube
 *     DRB - Background model cube
 *     DRG - Geometry factors cube
 *     DRX - Exposure map
 *
 * Each of the four files is mandatory.
 ***************************************************************************/
GCOMObservation::GCOMObservation(const GFilename& drename,
                                 const GFilename& drbname,
                                 const GFilename& drgname,
                                 const GFilename& drxname) : GObservation()
{
    // Initialise members
    init_members();

    // Load observation
    load(drename, drbname, drgname, drxname);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Unbinned observation constructor
 *
 * @param[in] evpname Event list FITS file name.
 * @param[in] timname Good Time Intervals FITS file name.
 * @param[in] oadnames List of Orbit Aspect Data FITS file names.
 *
 * Creates a COMPTEL unbinned observation by loading the event list, Good
 * Time Interval and Orbit Aspect Data from FITS files. All files are
 * mandatory.
 ***************************************************************************/
GCOMObservation::GCOMObservation(const GFilename&              evpname,
                                 const GFilename&              timname,
                                 const std::vector<GFilename>& oadnames)
{
    // Initialise members
    init_members();

    // Load observation
    load(evpname, timname, oadnames);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] obs COMPTEL observation.
 *
 * Creates COMPTEL observation by copying an existing COMPTEL observation.
 ***************************************************************************/
GCOMObservation::GCOMObservation(const GCOMObservation& obs) : GObservation(obs)
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
GCOMObservation::~GCOMObservation(void)
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
 * @param[in] obs COMPTEL observation.
 * @return COMPTEL observation.
 *
 * Assign COMPTEL observation to this object.
 ***************************************************************************/
GCOMObservation& GCOMObservation::operator=(const GCOMObservation& obs)
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
 * @brief Clear COMPTEL observation
 ***************************************************************************/
void GCOMObservation::clear(void)
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
 * @brief Clone COMPTEL observation
 *
 * @return Pointer to deep copy of COMPTEL observation.
 ***************************************************************************/
GCOMObservation* GCOMObservation::clone(void) const
{
    return new GCOMObservation(*this);
}


/***********************************************************************//**
 * @brief Set response function
 *
 * @param[in] rsp Response function.
 *
 * @exception GException::invalid_argument
 *            Specified response is not a COMPTEL response.
 *
 * Sets the response function for the observation.
 ***************************************************************************/
void GCOMObservation::response(const GResponse& rsp)
{
    // Get pointer on COM response
    const GCOMResponse* comrsp = dynamic_cast<const GCOMResponse*>(&rsp);
    if (comrsp == NULL) {
        std::string cls = std::string(typeid(&rsp).name());
        std::string msg = "Response of type \""+cls+"\" is "
                          "not a COMPTEL response. Please specify a COMPTEL "
                          "response as argument.";
        throw GException::invalid_argument(G_RESPONSE, msg);
    }

    // Clone response function
    m_response = *comrsp;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set response function
 *
 * @param[in] caldb Calibration database.
 * @param[in] rspname Name of COMPTEL response.
 *
 * Sets the response function by loading the response information from the
 * calibration database.
 ***************************************************************************/
void GCOMObservation::response(const GCaldb& caldb, const std::string& rspname)
{
    // Clear COM response function
    m_response.clear();

    // Set calibration database
    m_response.caldb(caldb);

    // Load instrument response function
    m_response.load(rspname);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read observation from XML element
 *
 * @param[in] xml XML element.
 *
 * Reads information for a COMPTEL observation from an XML element. The
 * method supports both an unbinned and a binned observation.
 *
 * For an unbinned observation the XML format is
 *
 *     <observation name="Crab" id="000001" instrument="COM">
 *       <parameter name="EVP" file="m16992_evp.fits"/>
 *       <parameter name="TIM" file="m10695_tim.fits"/>
 *       <parameter name="OAD" file="m20039_oad.fits"/>
 *       <parameter name="OAD" file="m20041_oad.fits"/>
 *       ...
 *     </observation>
 *
 * where the observation can contain an arbitrary number of OAD file
 * parameters. The @p file attribute provide either absolute or relative
 * file name. If a file name includes no access path it is assumed that
 * the file resides in the same location as the XML file.
 *
 * For a binned observation the XML format is
 *
 *     <observation name="Crab" id="000001" instrument="COM">
 *       <parameter name="DRE" file="m50438_dre.fits"/>
 *       <parameter name="DRB" file="m34997_drg.fits"/>
 *       <parameter name="DRG" file="m34997_drg.fits"/>
 *       <parameter name="DRX" file="m32171_drx.fits"/>
 *       <parameter name="IAQ" value="UNH(1.0-3.0)MeV"/>
 *     </observation>
 ***************************************************************************/
void GCOMObservation::read(const GXmlElement& xml)
{
    // Clear observation
    clear();

    // Extract instrument name
    m_instrument = xml.attribute("instrument");

    // If the XML elements has a "EVP" attribute we have an unbinned
    // observation
    if (gammalib::xml_has_par(xml, "EVP")) {

        // Get EVP and TIM file names
        std::string evpname = gammalib::xml_get_attr(G_READ, xml, "EVP", "file");
        std::string timname = gammalib::xml_get_attr(G_READ, xml, "TIM", "file");

        // Expand EVP and TIM file names
        evpname = gammalib::xml_file_expand(xml, evpname);
        timname = gammalib::xml_file_expand(xml, timname);

        // Get OAD file names
        std::vector<GFilename> oadnames;
        for (int i = 0; i < xml.elements("parameter"); ++i) {
            const GXmlElement* element = xml.element("parameter", i);
            if (element->attribute("name") == "OAD") {
                std::string oadname = element->attribute("file");
                oadname = gammalib::xml_file_expand(xml, oadname);
                oadnames.push_back(GFilename(oadname));
           }
        }

        // Load observation
        load(evpname, timname, oadnames);

    } // endif: unbinned observation

    // ... otherwise we have a binned observation
    else {

        // Get parameters
        std::string drename = gammalib::xml_get_attr(G_READ, xml, "DRE", "file");
        std::string drbname = gammalib::xml_get_attr(G_READ, xml, "DRB", "file");
        std::string drgname = gammalib::xml_get_attr(G_READ, xml, "DRG", "file");
        std::string drxname = gammalib::xml_get_attr(G_READ, xml, "DRX", "file");
        std::string iaqname = gammalib::xml_get_attr(G_READ, xml, "IAQ", "value");

        // Expand file names
        drename = gammalib::xml_file_expand(xml, drename);
        drbname = gammalib::xml_file_expand(xml, drbname);
        drgname = gammalib::xml_file_expand(xml, drgname);
        drxname = gammalib::xml_file_expand(xml, drxname);

        // Load observation
        load(drename, drbname, drgname, drxname);

        // Load IAQ
        response(GCaldb("cgro", "comptel"), iaqname);

    } // endelse: binned observation

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write observation into XML element
 *
 * @param[in] xml XML element.
 *
 * Writes information for a COMPTEL observation into an XML element. The
 * method supports both an unbinned and a binned observation.
 *
 * For an unbinned observation the XML format is
 *
 *     <observation name="Crab" id="000001" instrument="COM">
 *       <parameter name="EVP" file="m16992_evp.fits"/>
 *       <parameter name="TIM" file="m10695_tim.fits"/>
 *       <parameter name="OAD" file="m20039_oad.fits"/>
 *       <parameter name="OAD" file="m20041_oad.fits"/>
 *       ...
 *     </observation>
 *
 * where the observation can contain an arbitrary number of OAD file
 * parameters. The @p file attribute provide either absolute or relative
 * file name. If a file name includes no access path it is assumed that
 * the file resides in the same location as the XML file.
 *
 * For a binned observation the XML format is
 *
 *     <observation name="Crab" id="000001" instrument="COM">
 *       <parameter name="DRE" file="m50438_dre.fits"/>
 *       <parameter name="DRB" file="m34997_drg.fits"/>
 *       <parameter name="DRG" file="m34997_drg.fits"/>
 *       <parameter name="DRX" file="m32171_drx.fits"/>
 *       <parameter name="IAQ" value="UNH(1.0-3.0)MeV"/>
 *     </observation>
 ***************************************************************************/
void GCOMObservation::write(GXmlElement& xml) const
{
    // Allocate XML element pointer
    GXmlElement* par;

    // Handle unbinned observation
    if (is_unbinned()) {
    
        // Set EVP parameter
        par = gammalib::xml_need_par(G_WRITE, xml, "EVP");
        par->attribute("file", gammalib::xml_file_reduce(xml, m_evpname));

        // Set TIM parameter
        par = gammalib::xml_need_par(G_WRITE, xml, "TIM");
        par->attribute("file", gammalib::xml_file_reduce(xml, m_timname));

        // Remove all existing OAD parameters
        for (int i = 0; i < xml.elements("parameter"); ++i) {
            GXmlElement* element = xml.element("parameter", i);
            if (element->attribute("name") == "OAD") {
                xml.remove(i);
            }
        }

        // Set OAD parameters
        for (int i = 0; i < m_oadnames.size(); ++i) {
            par = static_cast<GXmlElement*>(xml.append(GXmlElement("parameter name=\"OAD\"")));
            par->attribute("file", gammalib::xml_file_reduce(xml, m_oadnames[i]));
        }

    } // endif: observation was unbinned

    // ... otherwise handle a binned observation
    else if (is_binned()) {

        // Set DRE parameter
        par = gammalib::xml_need_par(G_WRITE, xml, "DRE");
        par->attribute("file", gammalib::xml_file_reduce(xml, m_drename));

        // Set DRB parameter
        par = gammalib::xml_need_par(G_WRITE, xml, "DRB");
        par->attribute("file", gammalib::xml_file_reduce(xml, m_drbname));

        // Set DRG parameter
        par = gammalib::xml_need_par(G_WRITE, xml, "DRG");
        par->attribute("file", gammalib::xml_file_reduce(xml, m_drgname));

        // Set DRX parameter
        par = gammalib::xml_need_par(G_WRITE, xml, "DRX");
        par->attribute("file", gammalib::xml_file_reduce(xml, m_drxname));

        // Set IAQ parameter
        par = gammalib::xml_need_par(G_WRITE, xml, "DRX");
        par->attribute("value", m_response.rspname());

    } // endif: observation was binned

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load data for a binned observation
 *
 * @param[in] drename Event cube name.
 * @param[in] drbname Background cube name.
 * @param[in] drgname Geometry cube name.
 * @param[in] drxname Exposure map name.
 *
 * Load event cube from DRE file, background model from DRB file, geometry
 * factors from DRG file and the exposure map from the DRX file. All files
 * are mandatory.
 ***************************************************************************/
void GCOMObservation::load(const GFilename& drename,
                           const GFilename& drbname,
                           const GFilename& drgname,
                           const GFilename& drxname)
{
    // Load DRE
    load_dre(drename);

    // Load DRB
    load_drb(drbname);

    // Load DRG
    load_drg(drgname);

    // Load DRX
    load_drx(drxname);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load data for an unbinned observation
 *
 * @param[in] evpname Event list FITS file name.
 * @param[in] timname Good Time Intervals FITS file name.
 * @param[in] oadnames List of Orbit Aspect Data FITS file names.
 *
 * Loads the event list, Good Time Interval and Orbit Aspect Data for an
 * unbinned observation. All files are mandatory.
 ***************************************************************************/
void GCOMObservation::load(const GFilename&              evpname,
                           const GFilename&              timname,
                           const std::vector<GFilename>& oadnames)
{
    // Clear object
    clear();

    // Allocate event list
    GCOMEventList *evp = new GCOMEventList(evpname);
    m_events           = evp;

    // Extract observation information from event list
    // m_obs_id = ;
    // m_name   = ;

    // Extract time information from event list
    double tstart = m_events->gti().tstart().mjd();
    double tstop  = m_events->gti().tstop().mjd();
    m_ontime      = (tstop - tstart) * 86400.0;
    m_deadc       = 0.85;
    m_livetime    = m_deadc * m_ontime;

    // Load TIM data
    m_tim.load(timname);

    // Load OAD data
    for (int i = 0; i < oadnames.size(); ++i) {

        // Load OAD file
        GCOMOads oads(oadnames[i]);

        // Skip file if it is empty
        if (oads.size() < 1) {
            continue;
        }

        // Find index of where to insert OAD file
        int index = 0;
        for (; index < m_oads.size(); ++index) {
            if (oads[0].tstop() < m_oads[index][m_oads[index].size()-1].tstart()) {
                break;
            }
        }

        // Inserts Orbit Aspect Data
        if (index < m_oads.size()) {
            m_oads.insert(m_oads.begin()+index, oads);
        }
        else {
            m_oads.push_back(oads);
        }

    }

    // Store filenames
    m_evpname  = evpname;
    m_timname  = timname;
    m_oadnames = oadnames;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute event cube
 *
 * @param[in] dre COMPTEL event cube.
 *
 * @exception GException::invalid_argument
 *            Observation does not hold an event list.
 *            DRE cube has a non-positive Phibar bin size.
 *
 * Compute DRE event cube from event list (EVP), Good Time Intervals (TIM)
 * and Orbit Aspect Data (OAD).
 ***************************************************************************/
void GCOMObservation::compute_dre(GCOMDri& dre)
{
    // Debug
    #if defined(G_DEBUG_COMPUTE_DRE)
    std::cout << "GCOMObservation::compute_dre" << std::endl;
    std::cout << "============================" << std::endl;
    #endif

    // Get pointer to event list. Throw an exception if the observation
    // does not hold an event list
    const GCOMEventList* evp = dynamic_cast<const GCOMEventList*>(m_events);
    if (evp == NULL) {
        std::string msg = "Observation does not contain a COMPTEL event "
                          "list. Please load a COMPTEL event before calling "
                          "the method.";
        throw GException::invalid_argument(G_COMPUTE_DRE, msg);
    }

    // Check for positive Phibar bin size
    if (dre.phibin() <= 0.0) {
        std::string msg = "DRE cube has a non-positive Phibar bin size. Please "
                          "specify a DRE cube with a positive Phibar bin size.";
        throw GException::invalid_argument(G_COMPUTE_DRE, msg);
    }

    // Initialise variables
    int        last_tjd = 0;
    int        i_evt    = 0;
    GTime      tstart;
    GTime      tstop;
    GCOMStatus status;

    // Initialise statistics
    int num_used_superpackets    = 0;
    int num_skipped_superpackets = 0;
    int num_used_events          = 0;
    int num_event_outside_sp     = 0;
    int num_energy_too_low       = 0;
    int num_energy_too_high      = 0;
    int num_no_scatter_angle     = 0;
    int num_bad_modcom           = 0;
    int num_eha_too_small        = 0;
    int num_phibar_too_low       = 0;
    int num_phibar_too_high      = 0;
    int num_outside_dre          = 0;
    int num_outside_e1           = 0;
    int num_outside_e2           = 0;
    int num_outside_tof          = 0;
    int num_outside_psd          = 0;
    int num_bad_reflag           = 0;
    int num_event_before_dre     = 0;
    int num_event_after_dre      = 0;
    int num_d1module_off         = 0;
    int num_d2module_off         = 0;
    int num_processed            = 0;
    int num_superpackets         = 0;

    // Set all DRE bins to zero
    for (int i = 0; i < dre.size(); ++i) {
        dre[i] = 0.0;
    }

    // Signal that loop should be terminated
    bool terminate = false;

    // Loop over Orbit Aspect Data vector
    for (int i_oads = 0; i_oads < m_oads.size(); ++i_oads) {

        // Loop over Orbit Aspect Data
        for (int i_oad = 0; i_oad < m_oads[i_oads].size(); ++i_oad) {

            // Get reference to Orbit Aspect Data of superpacket
            const GCOMOad &oad = m_oads[i_oads][i_oad];

            // Increment superpacket counter
            num_superpackets++;

            // Skip superpacket if it is not within Good Time Intervals.
            // According to the original routine dalsel17.p7time.f only the
            // superpacket start time is checked.
            if (!m_tim.contains(oad.tstart())) {
                num_skipped_superpackets++;
                continue;
            }

            // Update superpackets statistics
            num_used_superpackets++;

            // Prepare Earth horizon angle comparison
            GSkyDir sky_geocentre;
            double  theta_geocentre = double(oad.gcel());
            double  phi_geocentre   = double(oad.gcaz());
            sky_geocentre.radec_deg(phi_geocentre, 90.0-theta_geocentre);

            // If TJD changed then update the module status table
            if (oad.tjd() != last_tjd) {
                last_tjd = oad.tjd();
                // TODO: Update module status table
            }

            // Update validity interval so that DRE will have the correct
            // interval
            if (num_used_superpackets == 1) {
                tstart = oad.tstart();
                tstop  = oad.tstop();
            }
            else {
                tstop  = oad.tstop();
            }

            // Collect all events within superpacket. Break if the end
            // of the event list was reached.
            for (; i_evt < evp->size(); ++i_evt) {

                // Get pointer to event
                const GCOMEventAtom* event = (*evp)[i_evt];

                // Break loop if the end of the superpacket was reached
                if (event->time() > oad.tstop()) {
                    break;
                }

                // Increase number of processed events
                num_processed++;

                // Skip event if it lies before the DRE start
                if (event->time() < dre.gti().tstart()) {
                    num_event_before_dre++;
                    continue;
                }

                // Break if event lies after the DRE stop
                if (event->time() > dre.gti().tstop()) {
                    num_event_after_dre++;
                    terminate = true;
                    break;
                }

                // Skip event if it lies before the superpacket start
                if (event->time() < oad.tstart()) {
                    num_event_outside_sp++;
                    continue;
                }

                // Skip event if it lies outside energy range
                if (event->energy() < dre.ebounds().emin()) {
                    num_energy_too_low++;
                    continue;
                }
                else if (event->energy() > dre.ebounds().emax()) {
                    num_energy_too_high++;
                    continue;
                }

                // Apply event selection, see PSSEVP.F
                // TODO: Hard coded for the moment, seems like EVP files
                //       contain already these thresholds !!!
                if (event->e1() < 0.070 || event->e1() > 20.0) {
                    num_outside_e1++;
                    continue;
                }
                if (event->e2() < 0.650 || event->e2() > 30.0) {
                    num_outside_e2++;
                    continue;
                }
                if (event->tof() < 115 || event->tof() > 130) {
                    num_outside_tof++;
                    continue;
                }
                if (event->psd() < 0 || event->psd() > 110) {
                    num_outside_psd++;
                    continue;
                }
                if (event->reflag() < 1) {
                    num_bad_reflag++;
                    continue;
                }

                // Check whether the event has a scatter angle determined.
                // This is signaled by a scatter angle -10e20 in radians.
                if (event->theta() < -1.0e3) {
                    num_no_scatter_angle++;
                    continue;
                }

                // Check for valid module IDs from MODCOM
                if (event->modcom() < 1 && event->modcom() > 98) {
                    num_bad_modcom++;
                    continue;
                }

                // Extract module IDs from MODCOM
                int id2 = (event->modcom()-1)/7 + 1;      // [1-14]
                int id1 =  event->modcom() - (id2-1) * 7; // [1-7]

                // Check module status. Skip all events where module is
                // signalled as not on.
                if (status.d1status(oad.tjd(), id1) != 1) {
                    num_d1module_off++;
                    continue;
                }
                if (status.d2status(oad.tjd(), id2) != 1) {
                    num_d2module_off++;
                    continue;
                }

                // Check minitelescope against that specified in DRI
                // definition and module status from database
                // TODO: IF ( BETAPQ(ID1,ID2)    .EQ. 1
                //           .AND.D1STAT(ID1) .NE. 0
                //           .AND.D2STAT(ID2) .NE. 0 ) THEN

                // Check the event is not falling in excluded region
                // of failed PMT
                // MPAR(10,11) = D2(x,y) of event in mm
                //
                // TODO: IF(D2STAT(ID2).EQ.1  .OR.
                // TODO  (D2STAT(ID2).EQ.9. .AND. EXD2R(ID2) .GE. 0.1. AND.
                // TODO (EXD2X(ID2)-MPAR(10)*.1)**2+(EXD2Y(ID2)-MPAR(11)*.1)**2
                // TODO .GT.      EXD2R(ID2)**2           )

                // Compute Compton scatter angle index. Skip if it's invalid.
                int iphibar = (event->phibar() - dre.phimin()) / dre.phibin();
                if (iphibar < 0) {
                    num_phibar_too_low++;
                    continue;
                }
                else if (iphibar >= dre.nphibar()) {
                    num_phibar_too_high++;
                    continue;
                }

                // Compute Earth horizon angle (do not trust the EVP value)
                GSkyDir sky_event;
                double  theta_event = double(event->theta());
                double  phi_event   = double(event->phi());
                sky_event.radec_deg(-phi_event, 90.0-theta_event);
                double eha = sky_event.dist_deg(sky_geocentre) - oad.georad();

                // Option: Earth horizon angle comparison
                #if defined(G_CHECK_EHA_COMPUTATION)
                if (std::abs(eha - event->eha()) > 1.5) {
                    std::string msg = "Earth horizon angle from EVP dataset ("+
                                      gammalib::str(event->eha())+" deg) "
                                      "differs from Earth horizon angle "
                                      "computed from Orbit Aspect Data ("+
                                      gammalib::str(eha)+" deg). Use the EVP "
                                      "value.";
                    gammalib::warning(G_COMPUTE_DRE, msg);
                }
                #endif

                // Check for Earth horizon angle. The limit is fixed here to
                // the values of 5, 7, ..., 53 deg for 25 Phibar layers that
                // was usually used for COMPTEL.
                // zeta = 5.0
                double ehamin = double(iphibar) * dre.phibin() + 5.0;
                if (event->eha() < ehamin) {
                    num_eha_too_small++;
                    continue;
                }

                // Now fill the DRE event array
                GSkyPixel pixel = dre.map().dir2pix(event->dir().dir());
                int       ichi  = int(pixel.x());
                int       ipsi  = int(pixel.y());
                if ((ichi >= 0) && (ichi < dre.nchi()) &&
                    (ipsi >= 0) && (ipsi < dre.npsi())) {
                    int inx   = ichi + (ipsi + iphibar * dre.npsi()) * dre.nchi();
                    dre[inx] += 1.0;
                    num_used_events++;
                }
                else {
                    num_outside_dre++;
                }

            } // endfor: collected events

            // Break if termination was signalled or if there are no more events
            if (terminate || i_evt >= evp->size()) {
                break;
            }

        } // endfor: looped over Orbit Aspect Data

        // Break if termination was signalled or if there are no more events
        if (terminate || i_evt >= evp->size()) {
            break;
        }

    } // endfor: looped over Orbit Aspect Data vector

    // Set Good Time interval for DRE
    dre.gti(GGti(tstart, tstop));

    // Debug
    #if defined(G_DEBUG_COMPUTE_DRE)
    std::cout << "Total number of superpackets .: " << num_superpackets << std::endl;
    std::cout << "Used superpackets ............: " << num_used_superpackets << std::endl;
    std::cout << "Skipped superpackets .........: " << num_skipped_superpackets << std::endl;
    std::cout << "Total number of events .......: " << evp->size() << std::endl;
    std::cout << "Processed events .............: " << num_processed << std::endl;
    std::cout << "Used events ..................: " << num_used_events << std::endl;
    std::cout << "Events outside superpacket ...: " << num_event_outside_sp << std::endl;
    std::cout << "Events before DRE GTI ........: " << num_event_before_dre << std::endl;
    std::cout << "Events after DRE GTI .........: " << num_event_after_dre << std::endl;
    std::cout << "Events with D1 module off ....: " << num_d1module_off << std::endl;
    std::cout << "Events with D2 module off ....: " << num_d2module_off << std::endl;
    std::cout << "Energy too low ...............: " << num_energy_too_low << std::endl;
    std::cout << "Energy too high ..............: " << num_energy_too_high << std::endl;
    std::cout << "Phibar too low ...............: " << num_phibar_too_low << std::endl;
    std::cout << "Phibar too high ..............: " << num_phibar_too_high << std::endl;
    std::cout << "Outside D1 selection .........: " << num_outside_e1 << std::endl;
    std::cout << "Outside D2 selection .........: " << num_outside_e2 << std::endl;
    std::cout << "Outside TOF selection ........: " << num_outside_tof << std::endl;
    std::cout << "Outside PSD selection ........: " << num_outside_psd << std::endl;
    std::cout << "Outside rejection flag sel. ..: " << num_bad_reflag << std::endl;
    std::cout << "Outside DRE cube .............: " << num_outside_dre << std::endl;
    std::cout << "No scatter angle .............: " << num_no_scatter_angle << std::endl;
    std::cout << "Bad module combination .......: " << num_bad_modcom << std::endl;
    std::cout << "Earth horizon angle too small : " << num_eha_too_small << std::endl;
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print observation information
 *
 * @param[in] chatter Chattiness.
 * @return String containing observation information.
 ***************************************************************************/
std::string GCOMObservation::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCOMObservation ===");

        // Append information
        result.append("\n"+gammalib::parformat("Name")+name());
        result.append("\n"+gammalib::parformat("Identifier")+id());
        result.append("\n"+gammalib::parformat("Instrument")+instrument());
        result.append("\n"+gammalib::parformat("Statistics")+statistics());
        result.append("\n"+gammalib::parformat("Ontime"));
        result.append(gammalib::str(ontime())+" sec");
        result.append("\n"+gammalib::parformat("Livetime"));
        result.append(gammalib::str(livetime())+" sec");
        result.append("\n"+gammalib::parformat("Deadtime correction"));
        result.append(gammalib::str(m_deadc));

        // Append response (if available)
        if (response()->rspname().length() > 0) {
            result.append("\n"+response()->print(gammalib::reduce(chatter)));
        }

        // Append events
        if (m_events != NULL) {
            result.append("\n"+m_events->print(gammalib::reduce(chatter)));
        }
        else {
            result.append("\n"+gammalib::parformat("Events")+"undefined");
        }

        // Append DRB, DRG and DRX
        //result.append("\n"+m_drb.print(chatter));
        //result.append("\n"+m_drg.print(chatter));
        //result.append("\n"+m_drx.print(chatter));

    } // endif: chatter was not silent

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                            Private methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GCOMObservation::init_members(void)
{
    // Initialise members
    m_instrument = "COM";
    m_response.clear();
    m_obs_id   = 0;
    m_ontime   = 0.0;
    m_livetime = 0.0;
    m_deadc    = 0.0;

    // Initialise members for binned observation
    m_drename.clear();
    m_drbname.clear();
    m_drgname.clear();
    m_drxname.clear();
    m_drb.clear();
    m_drg.clear();
    m_drx.clear();
    m_ewidth   = 0.0;

    // Initialise members for unbinned observation
    m_evpname.clear();
    m_timname.clear();
    m_oadnames.clear();
    m_tim.clear();
    m_oads.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] obs COMPTEL observation.
 ***************************************************************************/
void GCOMObservation::copy_members(const GCOMObservation& obs)
{
    // Copy members
    m_instrument = obs.m_instrument;
    m_response   = obs.m_response;
    m_obs_id     = obs.m_obs_id;
    m_ontime     = obs.m_ontime;
    m_livetime   = obs.m_livetime;
    m_deadc      = obs.m_deadc;

    // Copy members for binned observation
    m_drename    = obs.m_drename;
    m_drbname    = obs.m_drbname;
    m_drgname    = obs.m_drgname;
    m_drxname    = obs.m_drxname;
    m_drb        = obs.m_drb;
    m_drg        = obs.m_drg;
    m_drx        = obs.m_drx;
    m_ewidth     = obs.m_ewidth;

    // Copy members for unbinned observation
    m_evpname  = obs.m_evpname;
    m_timname  = obs.m_timname;
    m_oadnames = obs.m_oadnames;
    m_tim      = obs.m_tim;
    m_oads     = obs.m_oads;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCOMObservation::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Load event cube data from DRE file
 *
 * @param[in] drename DRE filename.
 *
 * Loads the event cube from a DRE file.
 ***************************************************************************/
void GCOMObservation::load_dre(const GFilename& drename)
{
    // Delete any existing event container (do not call clear() as we do not
    // want to delete the response function)
    if (m_events != NULL) delete m_events;
    m_events = NULL;

    // Allocate event cube
    m_events = new GCOMEventCube;

    // Open FITS file
    GFits fits(drename);

    // Read event cube
    m_events->read(fits);

    // Read observation attributes from primary extension
    GFitsHDU* hdu = fits[0];
    read_attributes(hdu);

    // Close FITS file
    fits.close();

    // Store event filename
    m_drename = drename;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load background model from DRB file
 *
 * @param[in] drbname DRB filename.
 *
 * @exception GException::invalid_value
 *            DRB data space incompatible with DRE data space.
 *
 * Load the background model from the primary image of the specified FITS
 * file.
 ***************************************************************************/
void GCOMObservation::load_drb(const GFilename& drbname)
{
    // Open FITS file
    GFits fits(drbname);

    // Get image
    const GFitsImage& image = *fits.image("Primary");

    // Load background model as sky map
    m_drb.read(image);

    // Correct WCS projection (HEASARC data format kluge)
    com_wcs_mer2car(m_drb);

    // Close FITS file
    fits.close();

    // Check map dimensions
    if (!check_map(m_drb)) {
        std::string msg = "DRB data cube \""+drbname+"\" incompatible with "
                          "DRE data cube \""+m_drename+"\".";
        throw GException::invalid_value(G_LOAD_DRB, msg);
    }

    // Store DRB filename
    m_drbname = drbname;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load geometry factors from DRG file
 *
 * @param[in] drgname DRG filename.
 *
 * @exception GException::invalid_value
 *            DRG data space incompatible with DRE data space.
 *
 * Load the geometry factors from the primary image of the specified FITS
 * file.
 ***************************************************************************/
void GCOMObservation::load_drg(const GFilename& drgname)
{
    // Open FITS file
    GFits fits(drgname);

    // Get image
    const GFitsImage& image = *fits.image("Primary");

    // Load geometry factors as sky map
    m_drg.read(image);

    // Correct WCS projection (HEASARC data format kluge)
    com_wcs_mer2car(m_drg);

    // Close FITS file
    fits.close();

    // Check map dimensions
    if (!check_map(m_drg)) {
        std::string msg = "DRG data cube \""+drgname+"\" incompatible with "
                          "DRE data cube \""+m_drename+"\".";
        throw GException::invalid_value(G_LOAD_DRB, msg);
    }

    // Store DRG filename
    m_drgname = drgname;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load exposure from DRX file
 *
 * @param[in] drxname DRX filename.
 *
 * Load the exposure map from the primary image of the specified FITS file.
 ***************************************************************************/
void GCOMObservation::load_drx(const GFilename& drxname)
{
    // Open FITS file
    GFits fits(drxname);

    // Get HDU
    const GFitsImage& image = *fits.image("Primary");

    // Load exposure map as sky map
    m_drx.read(image);

    // Correct WCS projection (HEASARC data format kluge)
    com_wcs_mer2car(m_drx);

    // Close FITS file
    fits.close();

    // Store DRX filename
    m_drxname = drxname;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Check if sky map is compatible with event cube
 *
 * @param[in] map Sky map.
 * @return True if map is compatible, false otherwise.
 *
 * Compares the dimension and the WCS definition of a sky map to that of the
 * event cube. If both are identical, true is returned, false otherwise.
 ***************************************************************************/
bool GCOMObservation::check_map(const GSkyMap& map) const
{
    // Get reference to event cube map
    const GSkyMap& ref = dynamic_cast<GCOMEventCube*>(m_events)->dre().map();

    // Compare dimensions
    bool same_dimension = ((map.nx()    == ref.nx()) &&
                           (map.ny()    == ref.ny()) &&
                           (map.nmaps() == ref.nmaps()));

    // Compare projections
    bool same_projection = (*(map.projection()) == *(ref.projection()));

    // Return
    return (same_dimension && same_projection);
}


/***********************************************************************//**
 * @brief Read observation attributes
 *
 * @param[in] hdu FITS HDU pointer
 *
 * Reads COM observation attributes from HDU. Mandatory attributes are
 *
 *     TSTART   - Start time (days)
 *     TSTOP    - Stop time (days)
 *     E_MIN    - Minimum energy (MeV)
 *     E_MAX    - Maximum energy (MeV)
 *
 * and optional attributes are
 *
 *     OBS_ID   - Observation identifier
 *     OBJECT   - Object
 *
 * Based on TSTART and TSTOP the ontime is computed. The deadtime fraction
 * is fixed to 15%, hence DEADC=0.85. The livetime is then computed by
 * multiplying the deadtime correction by the ontime, i.e.
 * LIVETIME = ONTIME * DEADC.
 *
 * Nothing is done if the HDU pointer is NULL.
 ***************************************************************************/
void GCOMObservation::read_attributes(const GFitsHDU* hdu)
{
    // Continue only if HDU is valid
    if (hdu != NULL) {

        // Read observation information
        m_obs_id = (hdu->has_card("OBS_ID")) ? hdu->real("OBS_ID") : 0;
        m_name   = (hdu->has_card("OBJECT")) ? hdu->string("OBJECT") : "unknown";

        // Compute ontime
        double tstart = hdu->real("TSTART");
        double tstop  = hdu->real("TSTOP");
        m_ontime   = (tstop - tstart) * 86400.0;
        m_deadc    = 0.85;
        m_livetime = m_deadc * m_ontime;

        // Compute energy width
        double emin = hdu->real("E_MIN");
        double emax = hdu->real("E_MAX");
        m_ewidth = emax - emin;

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
 *
 * @todo Implement method.
 ***************************************************************************/
void GCOMObservation::write_attributes(GFitsHDU* hdu) const
{
    // Continue only if HDU is valid
    if (hdu != NULL) {

    } // endif: HDU was valid

    // Return
    return;
}
