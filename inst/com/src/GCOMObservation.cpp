/***************************************************************************
 *             GCOMObservation.cpp - COMPTEL Observation class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2018 by Juergen Knoedlseder                         *
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
#include "GMath.hpp"
#include "GException.hpp"
#include "GFits.hpp"
#include "GCaldb.hpp"
#include "GSource.hpp"
#include "GModelSky.hpp"
#include "GModelSpectralConst.hpp"
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
#define G_LOAD_DRB                    "GCOMObservation::load_drb(GFilename&)"
#define G_LOAD_DRG                    "GCOMObservation::load_drg(GFilename&)"
#define G_ADD_DRM                        "GCOMObservation::add_drm(GSource&)"

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
 * @brief XML constructor
 *
 * @param[in] xml XML element.
 *
 * Constructs a COMPTEL observation from the information that is found in an
 * XML element.
 ***************************************************************************/
GCOMObservation::GCOMObservation(const GXmlElement& xml) : GObservation()
{
    // Initialise members
    init_members();

    // Read XML
    read(xml);

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
        par = gammalib::xml_need_par(G_WRITE, xml, "IAQ");
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

    // Initialise intermediate vector for OADs
    std::vector<GCOMOads> oads;

    // Load OAD data
    for (int i = 0; i < oadnames.size(); ++i) {

        // Load OAD file
        GCOMOads oad(oadnames[i]);

        // Skip file if it is empty
        if (oad.size() < 1) {
            continue;
        }

        // Find index of where to insert OAD file
        int index = 0;
        for (; index < oads.size(); ++index) {
            if (oad[0].tstop() < oads[index][oads[index].size()-1].tstart()) {
                break;
            }
        }

        // Inserts Orbit Aspect Data
        if (index < oads.size()) {
            oads.insert(oads.begin()+index, oad);
        }
        else {
            oads.push_back(oad);
        }

    }

    // Now put all OAD into a single container
    for (int i = 0; i < oads.size(); ++i) {
        m_oads.extend(oads[i]);
    }

    // Store filenames
    m_evpname  = evpname;
    m_timname  = timname;
    m_oadnames = oadnames;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return source model DRM cube
 *
 * @param[in] source Source.
 *
 * @return Source model DRM cube.
 ***************************************************************************/
const GCOMDri& GCOMObservation::drm(const GSource& source) const
{
    // Search source in DRM cache
    int index = 0;
    for (; index < m_drms.size(); ++index) {
        if (source.name() == m_drms[index].name()) {
            break;
        }
    }

    // If no source was found then add a cube to the DRM cache
    if (index >= m_drms.size()) {
        const_cast<GCOMObservation*>(this)->add_drm(source);
        index = m_drms.size() - 1;
    }

    // ... otherwise, if the spatial source model has changed then update
    // the DRM cache
    else {

        // Check if one of the source parameters has changed
        bool changed = (m_pars[index].size() != source.model()->size());
        if (!changed) {
            for (int i = 0; i < source.model()->size(); ++i) {
                if (m_pars[index][i] != (*(source.model()))[i].value()) {
                    changed = true;
                    break;
                }
            }
        }

        // If model has changed then recompute it
        if (changed) {

            // Setup sky model for DRM computation
            GModelSky model(*(source.model()), GModelSpectralConst());

            // Compute DRM
            const_cast<GCOMObservation*>(this)->m_drms[index].compute_drm((*this), model);

            // Store new model parameters in cache
            for (int i = 0; i < source.model()->size(); ++i) {
                const_cast<GCOMObservation*>(this)->m_pars[index][i] =
                    (*(source.model()))[i].value();
            }

        } // endif: source parameter(s) changed

    } // endelse: spatial source model has changed

    // Return reference to DRM cube
    return (m_drms[index]);
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
        result.append("\n"+gammalib::parformat("Statistic")+statistic());
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

    // Initialise members for response cache
    m_pars.clear();
    m_drms.clear();

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

    // Copy members for response cache
    m_pars = obs.m_pars;
    m_drms = obs.m_drms;

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
 *
 * The ontime is extracted from the Good Time Intervals. The deadtime
 * fraction is fixed to 15%, hence DEADC=0.85. The livetime is then computed
 * by multiplying the deadtime correction by the ontime, i.e.
 * LIVETIME = ONTIME * DEADC.
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

    // Extract ontime
    m_ontime = m_events->gti().ontime();

    // Set fixed deadtime fraction
    m_deadc = 0.85;

    // Compute livetime
    m_livetime = m_deadc * m_ontime;

    // Compute energy width
    m_ewidth = m_events->emax().MeV() - m_events->emin().MeV();

    // Read additional observation attributes from primary extension
    read_attributes(fits[0]);

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
 * @brief Add DRM cube to observation response cache
 *
 * @param[in] source Source.
 *
 * @exception GCOMObservation::add_drm
 *            No event cube specified.
 *
 * Adds a DRM cube based on the given source to the observation response
 * cache.
 ***************************************************************************/
void GCOMObservation::add_drm(const GSource& source)
{
    // Get pointer on COMPTEL event cube
    const GCOMEventCube* cube = dynamic_cast<const GCOMEventCube*>(m_events);
    if (cube == NULL) {
        std::string cls = std::string(typeid(m_events).name());
        std::string msg = "Events of type \""+cls+"\" is not a COMPTEL event "
                          "cube. Please specify a COMPTEL event cube when "
                          "using this method.";
        throw GException::invalid_argument(G_ADD_DRM, msg);
    }

    // Initialise DRM cube based on DRE cube
    GCOMDri drm = cube->dre();

    // Setup sky model for DRM computation
    GModelSky model(*(source.model()), GModelSpectralConst());

    // Compute DRM
    drm.compute_drm((*this), model);

    // Set model name
    drm.name(source.name());

    // Push model on stack
    m_drms.push_back(drm);

    // Push model parameters on stack
    std::vector<double> pars;
    for (int i = 0; i < source.model()->size(); ++i) {
        pars.push_back((*(source.model()))[i].value());
    }
    m_pars.push_back(pars);

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
 * Reads optional attributes are
 *
 *     OBS_ID   - Observation identifier
 *     OBJECT   - Object
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
