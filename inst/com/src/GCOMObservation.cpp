/***************************************************************************
 *             GCOMObservation.cpp - COMPTEL Observation class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2021 by Juergen Knoedlseder                         *
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
#include "GFitsHDU.hpp"
#include "GCaldb.hpp"
#include "GSource.hpp"
#include "GModels.hpp"
#include "GModelSky.hpp"
#include "GModelSpectralConst.hpp"
#include "GXmlElement.hpp"
#include "GCOMObservation.hpp"
#include "GCOMEventBin.hpp"
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
#define G_DRM                                "GCOMObservation::drm(GModels&)"
#define G_COMPUTE_DRB "GCOMObservation::compute_drb(std::string&, GCOMDri&, "\
                                                     "int&, int&, int&, int&"
#define G_COMPUTE_DRB_PHINOR  "GCOMObservation::compute_drb_phinor(GCOMDri&)"
#define G_COMPUTE_DRB_BGDLIXA         "GCOMObservation::compute_drb_bgdlixa("\
                                          "GCOMDri&, int&, int&, int&, int&)"
#define G_COMPUTE_DRB_BGDLIXE         "GCOMObservation::compute_drb_bgdlixe("\
                                          "GCOMDri&, int&, int&, int&, int&)"

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
 * @brief Binned observation DRI constructor
 *
 * @param[in] dre Event cube.
 * @param[in] drb Background cube.
 * @param[in] drg Geometry cube.
 * @param[in] drx Exposure map.
 *
 * Creates COMPTEL observation from DRI instances.
 ***************************************************************************/
GCOMObservation::GCOMObservation(const GCOMDri& dre,
                                 const GCOMDri& drb,
                                 const GCOMDri& drg,
                                 const GCOMDri& drx) : GObservation()
{
    // Initialise members
    init_members();

    // Store DRIs
    m_events = new GCOMEventCube(dre);
    m_drb    = drb;
    m_drg    = drg;
    m_drx    = drx;

    // Set attributes
    m_obs_id   = 0;
    m_ontime   = m_events->gti().ontime();
    m_deadc    = 0.965;
    m_livetime = m_deadc * m_ontime;
    m_ewidth   = m_events->emax().MeV() - m_events->emin().MeV();
    m_name     = "unknown";
    m_drename  = "";
    m_drbname  = "";
    m_drgname  = "";
    m_drxname  = "";

    // Return
    return;
}


/***********************************************************************//**
 * @brief Binned observation filename constructor
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

    // Extract observation information from event list
    GFits fits(evpname);
    read_attributes(fits[1]);
    fits.close();

    // Allocate event list
    GCOMEventList *evp = new GCOMEventList(evpname);
    m_events           = evp;

    // Load TIM data
    m_tim.load(timname);

    // Extract ontime from TIM and compute livetime assuming the value of
    // 0.965 from Rob van Dijk's thesis, page 62
    m_ontime   = m_tim.gti().ontime();
    m_deadc    = 0.965;
    m_livetime = m_deadc * m_ontime;

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
 * @brief Compute DRM cube
 *
 * @param[in] models Model container.
 * @return DRM cube (units of counts)
 *
 * @exception GException::invalid_value
 *            Observation does not contain an event cube.
 *
 * Computes a COMPTEL DRM cube from the information provided in a model
 * container. The values of the DRM cube are in units of counts.
 ***************************************************************************/
GCOMDri GCOMObservation::drm(const GModels& models) const
{
    // Get pointer on COMPTEL event cube
    const GCOMEventCube* cube = dynamic_cast<const GCOMEventCube*>(m_events);
    if (cube == NULL) {
        std::string cls = std::string(typeid(m_events).name());
        std::string msg = "Events of type \""+cls+"\" is not a COMPTEL event "
                          "cube. Please specify a COMPTEL event cube when "
                          "using this method.";
        throw GException::invalid_value(G_DRM, msg);
    }

    // Create DRM cube as copy of DRE cube
    GCOMEventCube drm_cube = *cube;

    // Get vector of model values
    GVector values = models.eval(*this);

    // Loop over all cube bins
    for (int i = 0; i < drm_cube.size(); ++i) {

        // Get pointer to cube bin
        GCOMEventBin* bin = drm_cube[i];

        // Compute model value for cube bin
        double model = values[i] * bin->size();

        // Store model value in cube bin
        bin->counts(model);

    } // endfor: looped over DRM cube bins

    // Get a copy of the DRM
    GCOMDri drm = drm_cube.dre();

    // Return DRM
    return drm;
}


/***********************************************************************//**
 * @brief Compute DRB cube
 *
 * @param[in] method Background method (PHINOR, BGDLIXA or BGDLIXE).
 * @param[in] drm DRM cube.
 * @param[in] nrunav BGDLIXA: number of bins used for running average.
 * @param[in] navgr BGDLIXA: number of bins used for averaging.
 * @param[in] nincl BGDLIXA: number of Phibar layers to include.
 * @param[in] nexcl BGDLIXA: number of Phibar layers to exclude.
 *
 * Computes a COMPTEL DRB cube using either the PHINOR or BGDLIXA method.
 * See the protected methods compute_drb_phinor() and compute_drb_bgdlixa()
 * for more information.
 ***************************************************************************/
void GCOMObservation::compute_drb(const std::string& method,
                                  const GCOMDri&     drm,
                                  const int&         nrunav,
                                  const int&         navgr,
                                  const int&         nincl,
                                  const int&         nexcl)
{
    // Branch to relevant method based on specified background method
    if (method == "PHINOR") {
        compute_drb_phinor(drm);
    }
    else if (method == "BGDLIXA") {
        compute_drb_bgdlixa(drm, nrunav, navgr, nincl, nexcl);
    }
    else if (method == "BGDLIXE") {
        compute_drb_bgdlixe(drm, nrunav, navgr, nincl, nexcl);
    }
    else {
        std::string msg = "Unknown background method \""+method+"\". "
                          "Specify either \"PHINOR\", \"BGDLIXA\" or "
                          "\"BGDLIXE\".";
        throw GException::invalid_argument(G_COMPUTE_DRB, msg);
    }

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

        // Append TIM (if available)
        if (m_tim.gti().size() > 0) {
            result.append("\n"+m_tim.print(gammalib::reduce(chatter)));
        }
        else {
            result.append("\n"+gammalib::parformat("TIM")+"undefined");
        }

        // Append OADs (if available)
        if (!m_oads.is_empty()) {
            result.append("\n"+m_oads.print(gammalib::reduce(chatter)));
        }
        else {
            result.append("\n"+gammalib::parformat("OADs")+"undefined");
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

    // Set fixed deadtime fraction (see Rob van Dijk's thesis, page 62)
    m_deadc = 0.965;

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

    // Read background model
    m_drb.read(image);

    // Correct WCS projection (HEASARC data format kluge)
    gammalib::com_wcs_mer2car(const_cast<GSkyMap&>(m_drb.map()));

    // Close FITS file
    fits.close();

    // Check map dimensions
    if (!check_dri(m_drb)) {
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

    // Read geometry factors
    m_drg.read(image);

    // Correct WCS projection (HEASARC data format kluge)
    gammalib::com_wcs_mer2car(const_cast<GSkyMap&>(m_drg.map()));

    // Close FITS file
    fits.close();

    // Check map dimensions
    if (!check_dri(m_drg)) {
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

    // Read exposure map
    m_drx.read(image);

    // Correct WCS projection (HEASARC data format kluge)
    gammalib::com_wcs_mer2car(const_cast<GSkyMap&>(m_drx.map()));

    // Close FITS file
    fits.close();

    // Store DRX filename
    m_drxname = drxname;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Check if DRI is compatible with event cube
 *
 * @param[in] dri DRI.
 * @return True if DRI is compatible, false otherwise.
 *
 * Compares the dimension and the WCS definition of a DRI to that of the
 * event cube. If both are identical, true is returned, false otherwise.
 ***************************************************************************/
bool GCOMObservation::check_dri(const GCOMDri& dri) const
{
    // Get references to event cube map and DRI map
    const GSkyMap& ref = dynamic_cast<GCOMEventCube*>(m_events)->dre().map();
    const GSkyMap& map = dri.map();

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


/***********************************************************************//**
 * @brief Compute DRB cube using PHINOR method
 *
 * @param[in] drm DRM cube.
 *
 * @exception GException::invalid_value
 *            Observation does not contain an event cube
 * @exception GException::invalid_argument
 *            DRM cube is incompatible with DRE
 ***************************************************************************/
void GCOMObservation::compute_drb_phinor(const GCOMDri& drm)
{
    // Extract COMPTEL event cube
    const GCOMEventCube* cube = dynamic_cast<const GCOMEventCube*>(this->events());
    if (cube == NULL) {
        std::string msg = "Observation \""+this->name()+"\" ("+this->id()+") "
                          "does not contain a COMPTEL event cube. Please "
                          "specify a COMPTEL observation containing and event "
                          "cube.";
        throw GException::invalid_value(G_COMPUTE_DRB_PHINOR, msg);
    }

    // Check DRM
    if (!check_dri(drm)) {
        std::string msg = "Specified DRM cube is incompatible with DRE. Please "
                          "specify a DRM with a data-space definition that is "
                          "identical to that of the DRE.";
        throw GException::invalid_argument(G_COMPUTE_DRB_PHINOR, msg);
    }

    // Initialise DRB by cloning DRG
    m_drb = m_drg;

    // Get DRI sky maps
    const GSkyMap& map_dre = cube->dre().map();
    const GSkyMap& map_drm = drm.map();
    GSkyMap        map_drg = get_weighted_drg_map();
    GSkyMap&       map_drb = const_cast<GSkyMap&>(m_drb.map());

    // Get data space dimensions
    int npix    = m_drg.map().npix();
    int nphibar = m_drg.nphibar();

    // Phibar normalise DRB
    for (int iphibar = 0; iphibar < nphibar; ++iphibar) {
        double sum_dre = 0.0;
        double sum_drm = 0.0;
        double sum_drg = 0.0;
        for (int ipix = 0; ipix < npix; ++ipix) {
            sum_dre += map_dre(ipix, iphibar);
            sum_drm += map_drm(ipix, iphibar);
            sum_drg += map_drg(ipix, iphibar);
        }
        if (sum_drg > 0.0) {
            double norm = (sum_dre - sum_drm) / sum_drg;
            for (int ipix = 0; ipix < npix; ++ipix) {
                map_drb(ipix, iphibar) = map_drg(ipix, iphibar) * norm;
            }
        }
        else {
            for (int ipix = 0; ipix < npix; ++ipix) {
                map_drb(ipix, iphibar) = 0.0;
            }
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute DRB cube using BGDLIXA method
 *
 * @param[in] drm DRM cube.
 * @param[in] nrunav Number of bins used for running average.
 * @param[in] navgr Number of bins used for averaging.
 * @param[in] nincl Number of Phibar layers to include.
 * @param[in] nexcl Number of Phibar layers to exclude.
 *
 * @exception GException::invalid_value
 *            Observation does not contain an event cube
 * @exception GException::invalid_argument
 *            DRM cube is incompatible with DRE
 *
 * Computes a DRB cube using the BGDLIXA method that is documented in Rob
 * van Dijk's PhD thesis. The revelant equations from the thesis that are
 * implemented here are 3.12, 3.12 and 3.14.
 ***************************************************************************/
void GCOMObservation::compute_drb_bgdlixa(const GCOMDri& drm,
                                          const int&     nrunav,
                                          const int&     navgr,
                                          const int&     nincl,
                                          const int&     nexcl)
{
    // Extract COMPTEL event cube
    const GCOMEventCube* cube = dynamic_cast<const GCOMEventCube*>(this->events());
    if (cube == NULL) {
        std::string msg = "Observation \""+this->name()+"\" ("+this->id()+") "
                          "does not contain a COMPTEL event cube. Please "
                          "specify a COMPTEL observation containing and event "
                          "cube.";
        throw GException::invalid_value(G_COMPUTE_DRB_BGDLIXA, msg);
    }

    // Check DRM
    if (!check_dri(drm)) {
        std::string msg = "Specified DRM cube is incompatible with DRE. Please "
                          "specify a DRM with a data-space definition that is "
                          "identical to that of the DRE.";
        throw GException::invalid_argument(G_COMPUTE_DRB_BGDLIXA, msg);
    }

    // Initialise DRB and scratch DRI by cloning DRG
    m_drb       = m_drg;
    GCOMDri dri = m_drg;

    // Get references to relevant sky maps
    const GSkyMap& map_dre = cube->dre().map();
    const GSkyMap& map_drm = drm.map();
    GSkyMap        map_drg = get_weighted_drg_map();
    GSkyMap&       map_drb = const_cast<GSkyMap&>(m_drb.map());
    GSkyMap&       map_dri = const_cast<GSkyMap&>(dri.map());

    // Get data space dimensions
    int nchi    = m_drg.nchi();
    int npsi    = m_drg.npsi();
    int nphibar = m_drg.nphibar();
    int npix    = nchi * npsi;
    int nsize   = npix * nphibar;

    // Precompute half running average lengths
    int navgr2 = int(navgr/2);
    int nexcl2 = int(nexcl/2);
    int nincl2 = int(nincl/2);

    // Initialise DRB and scratch DRI
    for (int iphibar = 0; iphibar < nphibar; ++iphibar) {
        for (int ipix = 0; ipix < npix; ++ipix) {
            map_drb(ipix, iphibar) = 0.0;
            map_dri(ipix, iphibar) = 0.0;
        }
    }

    // Phibar normalise DRB (Equation 3.12 in Rob van Dijk's PhD thesis).
    for (int iphibar = 0; iphibar < nphibar; ++iphibar) {
        double sum_dre = 0.0;
        double sum_drm = 0.0;
        double sum_drg = 0.0;
        for (int ipix = 0; ipix < npix; ++ipix) {
            sum_dre += map_dre(ipix, iphibar);
            sum_drm += map_drm(ipix, iphibar);
            sum_drg += map_drg(ipix, iphibar);
        }
        if (sum_drg > 0.0) {
            double norm = (sum_dre - sum_drm) / sum_drg;
            for (int ipix = 0; ipix < npix; ++ipix) {
                map_dri(ipix, iphibar) = map_drg(ipix, iphibar) * norm;
            }
        }
    }

    // Do 3D running average of DRE and Phibar normalised DRG (Equation 3.13
    // in Rob van Dijk's PhD thesis). The result is a Phibar normalised DRG
    // that is normalised over the 3D region to the DRE. The 3D region is all
    // Phibar layers and [-nrunav,+nrunav] Chi/Psi pixels around the pixel of
    // consideration.
    if (nrunav >= 1) {

        // Loop over Chi pixels
        for (int ichi = 0; ichi < nchi; ++ichi) {

            // Compute running average window in Chi
            int kchi_min = ichi - nrunav;
            int kchi_max = ichi + nrunav;
            if (kchi_min < 0) {
                kchi_min = 0;
            }
            if (kchi_max >= nchi) {
                kchi_max = nchi - 1;
            }

            // Loop over Psi pixels
            for (int ipsi = 0; ipsi < npsi; ++ipsi) {

                // Compute running average window in Psi
                int kpsi_min = ipsi - nrunav;
                int kpsi_max = ipsi + nrunav;
                if (kpsi_min < 0) {
                    kpsi_min = 0;
                }
                if (kpsi_max >= npsi) {
                    kpsi_max = npsi - 1;
                }

                // Initialise sums
                double sum_dre = 0.0;
                double sum_drm = 0.0;
                double sum_dri = 0.0;

                // Compute running average
                for (int kchi = kchi_min; kchi <= kchi_max; ++kchi) {
                    for (int kpsi = kpsi_min; kpsi <= kpsi_max; ++kpsi) {
                        int kpix = kchi + kpsi * nchi;
                        for (int kphibar = 0; kphibar < nphibar; ++kphibar) {
                            if (map_drg(kpix, kphibar) != 0.0) {
                                sum_dre += map_dre(kpix, kphibar);
                                sum_drm += map_drm(kpix, kphibar);
                                sum_dri += map_dri(kpix, kphibar);
                            }
                        }
                    }
                }

                // Renormalise scratch array
                if (sum_dri != 0.0) {
                    int    ipix = ichi + ipsi * nchi;
                    double norm = (sum_dre - sum_drm) / sum_dri;
                    for (int iphibar = 0; iphibar < nphibar; ++iphibar) {
                        map_dri(ipix, iphibar) *= norm;
                    }
                }

            } // endfor: looped over Phi

        } // endfor: looped over Chi

    } // endif: running averaging was requested

    // First part of Equation (3.14) in Rob van Dijk's PhD thesis.
    // (DRB = B^0C_L / sum B^0C_L). The DRB is pre-computed by adjusting the
    // scratch DRI over nincl Phibar layers.
    int ksel1 = 0;
    int kex1  = 0;
    int kex2  = 0;
    int ksel2 = 0;
    for (int iphibar = 0; iphibar < nphibar; ++iphibar) {

        // Get Phibar index range for sums
        get_bgdlixa_phibar_indices(iphibar, nincl, nexcl,
                                   &ksel1, &kex1, &kex2, &ksel2);

        // Loop over Chi/Psi pixels
        for (int ipix = 0; ipix < npix; ++ipix) {

            // Initialise DRB value
            map_drb(ipix, iphibar) = 0.0;

            // Continue only if scratch DRI is non-zero
            if (map_dri(ipix, iphibar) != 0.0) {

                // Initialise sum
                double sum_dri = 0.0;

                // Take sum over [ksel1, kex1]
                if (kex1 >= 0) {
                    for (int kphibar = ksel1; kphibar <= kex1; ++kphibar) {
                        sum_dri += map_dri(ipix, kphibar);
                    }
                }

                // Take sum over [kxe2, ksel2]
                if (kex2 < nphibar) {
                    for (int kphibar = kex2; kphibar <= ksel2; ++kphibar) {
                        sum_dri += map_dri(ipix, kphibar);
                    }
                }

                // If sum is not zero then set DRB
                if (sum_dri != 0.0) {
                    map_drb(ipix, iphibar) = map_dri(ipix, iphibar) / sum_dri;
                }

            } // endif: scratch DRI was non-zero

        } // endfor: looped over Chi/Psi pixels

    } // endfor: looped over Phibar

    // Second part of Equation (3.14) in Rob van Dijk's PhD thesis that
    // corresponds to G * [E / G]^s that will be stored in scratch DRI.
    //
    // NOTES:
    // - we replaced DRE by DRE-DRM in the nominator of the running
    //   average since and source component should of course be
    //   subtracted
    // - we implemented an alternative running average computation where
    //   the ratio between the DRE and DRG sums is taken. This avoids
    //   the divergence of the ratio at the edge of the DRG.
    for (int ichi = 0; ichi < nchi; ++ichi) {

        // Compute running average window in Chi
        int kchi_min = ichi - navgr2;
        int kchi_max = ichi + navgr2;
        if (kchi_min < 0) {
            kchi_min = 0;
        }
        if (kchi_max >= nchi) {
            kchi_max = nchi - 1;
        }

        // Loop over Psi pixels
        for (int ipsi = 0; ipsi < npsi; ++ipsi) {

            // Compute running average window in Psi
            int kpsi_min = ipsi - navgr2;
            int kpsi_max = ipsi + navgr2;
            if (kpsi_min < 0) {
                kpsi_min = 0;
            }
            if (kpsi_max >= npsi) {
                kpsi_max = npsi - 1;
            }

            // Loop over Phibar layers
            for (int iphibar = 0; iphibar < nphibar; ++iphibar) {

                // Initialise sums
                double sum_dre = 0.0;
                double sum_drm = 0.0;
                double sum_drg = 0.0;

                // Compute running average
                for (int kchi = kchi_min; kchi <= kchi_max; ++kchi) {
                    for (int kpsi = kpsi_min; kpsi <= kpsi_max; ++kpsi) {
                        int kpix = kchi + kpsi * nchi;
                        sum_dre += map_dre(kpix, iphibar);
                        sum_drm += map_drm(kpix, iphibar);
                        sum_drg += map_drg(kpix, iphibar);
                    }
                }

                // Multiply average by DRG and store it in scratch DRI
                int ipix = ichi + ipsi * nchi;
                if (sum_drg > 0.0) {
                    double norm = (sum_dre - sum_drm) / sum_drg;
                    map_dri(ipix, iphibar) = map_drg(ipix, iphibar) * norm;
                }
                else {
                    map_dri(ipix, iphibar) = 0.0;
                }

            } // endfor: looped over Phibar layers

        } // endfor: looped over Psi pixels

    } // endfor: looped over Chi pixels

    // Third part of Equation (3.14) in Rob van Dijk's PhD thesis that sums
    // G * [E / G]^s over a subset of Phibar layers (the first parenthesis)
    // and multplies with B^0C_L / sum B^0C_L
    for (int iphibar = 0; iphibar < nphibar; ++iphibar) {

        // Get Phibar index range for sums
        get_bgdlixa_phibar_indices(iphibar, nincl, nexcl,
                                   &ksel1, &kex1, &kex2, &ksel2);

        // Loop over Chi/Psi pixels
        for (int ipix = 0; ipix < npix; ++ipix) {

            // If DRB is not empty then correct it
            if (map_drb(ipix, iphibar) != 0.0) {

                // Initialise DRI sum
                double sum_dri = 0.0;

                // Take sum over [ksel1, kex1]
                if (kex1 >= 0) {
                    for (int kphibar = ksel1; kphibar <= kex1; ++kphibar) {
                        sum_dri += map_dri(ipix, kphibar);
                    }
                }

                // Take sum over [kxe2, ksel2]
                if (kex2 < nphibar) {
                    for (int kphibar = kex2; kphibar <= ksel2; ++kphibar) {
                        sum_dri += map_dri(ipix, kphibar);
                    }
                }

                // Correct DRB
                map_drb(ipix, iphibar) *= sum_dri;

            } // endif: DRB was not empty

        } // endfor: looped over Chi/Psi pixels

    } // endfor: looped over Phibar

    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute DRB cube using BGDLIXE method
 *
 * @param[in] drm DRM cube.
 * @param[in] nrunav Number of bins used for running average.
 * @param[in] navgr Number of bins used for averaging.
 * @param[in] nincl Number of Phibar layers to include.
 * @param[in] nexcl Number of Phibar layers to exclude.
 *
 * @exception GException::invalid_value
 *            Observation does not contain an event cube
 * @exception GException::invalid_argument
 *            DRM cube is incompatible with DRE
 *
 * Computes a DRB cube using the BGDLIXE method. This method differs from
 * the BGDLIXA method in the last step.
 ***************************************************************************/
void GCOMObservation::compute_drb_bgdlixe(const GCOMDri& drm,
                                          const int&     nrunav,
                                          const int&     navgr,
                                          const int&     nincl,
                                          const int&     nexcl)
{
    // Extract COMPTEL event cube
    const GCOMEventCube* cube = dynamic_cast<const GCOMEventCube*>(this->events());
    if (cube == NULL) {
        std::string msg = "Observation \""+this->name()+"\" ("+this->id()+") "
                          "does not contain a COMPTEL event cube. Please "
                          "specify a COMPTEL observation containing and event "
                          "cube.";
        throw GException::invalid_value(G_COMPUTE_DRB_BGDLIXE, msg);
    }

    // Check DRM
    if (!check_dri(drm)) {
        std::string msg = "Specified DRM cube is incompatible with DRE. Please "
                          "specify a DRM with a data-space definition that is "
                          "identical to that of the DRE.";
        throw GException::invalid_argument(G_COMPUTE_DRB_BGDLIXE, msg);
    }

    // Initialise DRB and scratch DRI by cloning DRG
    m_drb       = m_drg;
    GCOMDri dri = m_drg;

    // Get references to relevant sky maps
    const GSkyMap& map_dre = cube->dre().map();
    const GSkyMap& map_drm = drm.map();
    GSkyMap        map_drg = get_weighted_drg_map();
    GSkyMap&       map_drb = const_cast<GSkyMap&>(m_drb.map());
    GSkyMap&       map_dri = const_cast<GSkyMap&>(dri.map());

    // Get data space dimensions
    int nchi    = m_drg.nchi();
    int npsi    = m_drg.npsi();
    int nphibar = m_drg.nphibar();
    int npix    = nchi * npsi;
    int nsize   = npix * nphibar;

    // Precompute half running average lengths
    int navgr2 = int(navgr/2);
    int nexcl2 = int(nexcl/2);
    int nincl2 = int(nincl/2);

    // Initialise DRB and scratch DRI
    for (int iphibar = 0; iphibar < nphibar; ++iphibar) {
        for (int ipix = 0; ipix < npix; ++ipix) {
            map_drb(ipix, iphibar) = 0.0;
            map_dri(ipix, iphibar) = 0.0;
        }
    }

    // Phibar normalise DRG (Equation 3.12 in Rob van Dijk's PhD thesis).
    for (int iphibar = 0; iphibar < nphibar; ++iphibar) {
        double sum_dre = 0.0;
        double sum_drm = 0.0;
        double sum_drg = 0.0;
        for (int ipix = 0; ipix < npix; ++ipix) {
            sum_dre += map_dre(ipix, iphibar);
            sum_drm += map_drm(ipix, iphibar);
            sum_drg += map_drg(ipix, iphibar);
        }
        if (sum_drg > 0.0) {
            double norm = (sum_dre - sum_drm) / sum_drg;
            for (int ipix = 0; ipix < npix; ++ipix) {
                map_dri(ipix, iphibar) = map_drg(ipix, iphibar) * norm;
            }
        }
    }

    // Fit Phibar normalised DRG to DRE-DRM
    for (int ichi = 0; ichi < nchi; ++ichi) {

        // Compute running average window in Chi
        int kchi_min = ichi - navgr2;
        int kchi_max = ichi + navgr2;
        if (kchi_min < 0) {
            kchi_min = 0;
        }
        if (kchi_max >= nchi) {
            kchi_max = nchi - 1;
        }

        // Loop over Psi pixels
        for (int ipsi = 0; ipsi < npsi; ++ipsi) {

            // Compute running average window in Psi
            int kpsi_min = ipsi - navgr2;
            int kpsi_max = ipsi + navgr2;
            if (kpsi_min < 0) {
                kpsi_min = 0;
            }
            if (kpsi_max >= npsi) {
                kpsi_max = npsi - 1;
            }

            // Loop over Phibar layers
            for (int iphibar = 0; iphibar < nphibar; ++iphibar) {

                 // Get Phibar index range for sums
                int ksel1 = 0;
                int kex1  = 0;
                int kex2  = 0;
                int ksel2 = 0;
                get_bgdlixa_phibar_indices(iphibar, nincl, nexcl,
                                           &ksel1, &kex1, &kex2, &ksel2);

                // Initialise sums
                double sum_dre = 0.0;
                double sum_drm = 0.0;
                double sum_dri = 0.0;

                // Compute running average
                for (int kchi = kchi_min; kchi <= kchi_max; ++kchi) {
                    for (int kpsi = kpsi_min; kpsi <= kpsi_max; ++kpsi) {

                        // Get index
                        int kpix = kchi + kpsi * nchi;

                        // Take sum over [ksel1, kex1]
                        if (kex1 >= 0) {
                            for (int kphibar = ksel1; kphibar <= kex1; ++kphibar) {
                                sum_dre += map_dre(kpix, kphibar);
                                sum_drm += map_drm(kpix, kphibar);
                                sum_dri += map_dri(kpix, kphibar);
                            }
                        }

                        // Take sum over [kxe2, ksel2]
                        if (kex2 < nphibar) {
                            for (int kphibar = kex2; kphibar <= ksel2; ++kphibar) {
                                sum_dre += map_dre(kpix, kphibar);
                                sum_drm += map_drm(kpix, kphibar);
                                sum_dri += map_dri(kpix, kphibar);
                            }
                        }

                    } // endfor: looped over kpsi
                } // endfor: looped over kchi

                // Renormalize DRI locally to DRE-DRM
                int ipix = ichi + ipsi * nchi;
                if (sum_dri > 0.0) {
                    double norm            = (sum_dre - sum_drm) / sum_dri;
                    map_drb(ipix, iphibar) = map_dri(ipix, iphibar) * norm;
                }
                else {
                    map_drb(ipix, iphibar) = 0.0;
                }

            } // endfor: looped over Phibar layers

        } // endfor: looped over Psi pixels

    } // endfor: looped over Chi pixels

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return weighted DRG map
 *
 * @return Weighted DRG map.
 *
 * Returns a DRG as sky map where each pixel of the map is multiplied by the
 * solidangle of the pixel.
 ***************************************************************************/
GSkyMap GCOMObservation::get_weighted_drg_map(void) const
{
    // Initialise
    GSkyMap drg = m_drg.map();

    // Get data space dimensions
    int npix    = drg.npix();
    int nphibar = drg.nmaps();

    // Weight map
    for (int ipix = 0; ipix < npix; ++ipix) {
        double omega = drg.solidangle(ipix);
        for (int iphibar = 0; iphibar < nphibar; ++iphibar) {
            drg(ipix,iphibar) *= omega;
        }
    }

    // Return weighted DRG map
    return drg;
}


/***********************************************************************//**
 * @brief Compute Phibar index range for BGDLIXA background method
 *
 * @param[in] iphibar Phibar layer index.
 * @param[in] nincl Number of Phibar layers to include.
 * @param[in] nexcl Number of Phibar layers to exclude.
 * @param[out] isel1 Start index for first sum.
 * @param[out] iex1 Stop index for first sum.
 * @param[out] iex2 Start index for second sum.
 * @param[out] isel2 Stop index for second sum.
 *
 * This method is a helper method for the BGDLIXA background method. It
 * computes the Phibar index range for the third step of the background
 * model computation. The third step corresponds to Equation (3.14) in Rob
 * van Dijk's PhD thesis.
 *
 * The Phibar sum will be taken over the index ranges [isel1, iex1] and
 * [ixe2, isel2].
 ***************************************************************************/
void GCOMObservation::get_bgdlixa_phibar_indices(const int& iphibar,
                                                 const int& nincl,
                                                 const int& nexcl,
                                                 int*       isel1,
                                                 int*       iex1,
                                                 int*       iex2,
                                                 int*       isel2) const
{
    // Get number of Phibar layers
    int nphibar = m_drg.nphibar();
    int imax    = nphibar - 1;

    // Precompute half running average lengths
    int nexcl2 = int(nexcl/2);
    int nincl2 = int(nincl/2);

    // Compute Phibar index range without respecting boundaries
    *iex1  = iphibar - nexcl2 - 1;
    *iex2  = iphibar + nexcl2 + 1;
    *isel1 = iphibar - nincl2;
    *isel2 = iphibar + nincl2;

    // If no Phibar layers are to be excluded then make sure that the
    // Phibar layer range is continuous
    if (nexcl == 0) {
        (*iex2)--;
    }

    // Make sure that start index of first sum and stop index of second
    // some are within the validity range
    if (*isel1 < 0) {
        *isel1 = 0;
    }
    if (*isel2 > imax) {
        *isel2 = imax;
    }
 
    // Make sure that the sums use always ncincl Phibar layers, also near
    // the edges
    if (*isel1 == 0) {
        *isel2 = nincl - 1;
        if (nincl > imax) {
            *isel2 = imax;
        }

    }
    if (*isel2 == imax) {
        *isel1 = nphibar - nincl;
        if (*isel1 < 0) {
            *isel1 = 0;
        }
    }

    // Return
    return;
}
