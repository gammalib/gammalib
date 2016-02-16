/***************************************************************************
 *             GCOMObservation.cpp - COMPTEL Observation class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2016 by Juergen Knoedlseder                         *
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
 * @brief Instrument constructor
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
 * @param[in] rspname Name of COMPTEL response.
 * @param[in] caldb Calibration database.
 *
 * Sets the response function by loading the response information from the
 * calibration database.
 ***************************************************************************/
void GCOMObservation::response(const std::string& rspname, const GCaldb& caldb)
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
 * Reads information for a COMPTEL observation from an XML element.
 * The expected format of the XML element is
 *
 *     <observation name="..." id="..." instrument="...">
 *       <parameter name="DRE" file="..."/>
 *       <parameter name="DRB" file="..."/>
 *       <parameter name="DRG" file="..."/>
 *       <parameter name="DRX" file="..."/>
 *       <parameter name="IAQ" value="..."/>
 *     </observation>
 *
 * for a binned observation.
 ***************************************************************************/
void GCOMObservation::read(const GXmlElement& xml)
{
    // Clear observation
    clear();

    // Extract instrument name
    m_instrument = xml.attribute("instrument");

    // Get parameters
    std::string drename = gammalib::xml_get_attr(G_READ, xml, "DRE", "file");
    std::string drbname = gammalib::xml_get_attr(G_READ, xml, "DRB", "file");
    std::string drgname = gammalib::xml_get_attr(G_READ, xml, "DRG", "file");
    std::string drxname = gammalib::xml_get_attr(G_READ, xml, "DRX", "file");
    std::string iaqname = gammalib::xml_get_attr(G_READ, xml, "IAQ", "value");

    // Load observation
    load(drename, drbname, drgname, drxname);

    // Load IAQ
    response(iaqname, GCaldb("cgro", "comptel"));

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write observation into XML element
 *
 * @param[in] xml XML element.
 *
 * Writes information for a COMPTEL observation into an XML element. The
 * expected format of the XML element is
 *
 *     <observation name="..." id="..." instrument="...">
 *       <parameter name="DRE" file="..."/>
 *       <parameter name="DRB" file="..."/>
 *       <parameter name="DRG" file="..."/>
 *       <parameter name="DRX" file="..."/>
 *       <parameter name="IAQ" value="..."/>
 *     </observation>
 *
 * for a binned observation.
 ***************************************************************************/
void GCOMObservation::write(GXmlElement& xml) const
{
    // Allocate XML element pointer
    GXmlElement* par;

    // Set DRE parameter
    par = gammalib::xml_need_par(G_WRITE, xml, "DRE");
    par->attribute("file", m_drename);

    // Set DRB parameter
    par = gammalib::xml_need_par(G_WRITE, xml, "DRB");
    par->attribute("file", m_drbname);

    // Set DRG parameter
    par = gammalib::xml_need_par(G_WRITE, xml, "DRG");
    par->attribute("file", m_drgname);

    // Set DRX parameter
    par = gammalib::xml_need_par(G_WRITE, xml, "DRX");
    par->attribute("file", m_drxname);

    // Set IAQ parameter
    par = gammalib::xml_need_par(G_WRITE, xml, "DRX");
    par->attribute("value", m_response.rspname());

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load data for binned analysis
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
 * @brief Print observation information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
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
        result.append("\n"+gammalib::parformat("Energy band"));
        result.append(gammalib::str(ewidth())+" MeV");

        // Append pointing
        result.append("\n"+m_pointing.print(gammalib::reduce(chatter)));

        // Append response
        result.append("\n"+response()->print(gammalib::reduce(chatter)));

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
    m_drename.clear();
    m_drbname.clear();
    m_drgname.clear();
    m_drxname.clear();
    m_drb.clear();
    m_drg.clear();
    m_drx.clear();
    m_pointing.clear();
    m_response.clear();
    m_obs_id   = 0;
    m_ontime   = 0.0;
    m_livetime = 0.0;
    m_deadc    = 0.0;
    m_ewidth   = 0.0;

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
    m_drename    = obs.m_drename;
    m_drbname    = obs.m_drbname;
    m_drgname    = obs.m_drgname;
    m_drxname    = obs.m_drxname;
    m_drb        = obs.m_drb;
    m_drg        = obs.m_drg;
    m_drx        = obs.m_drx;
    m_response   = obs.m_response;
    m_pointing   = obs.m_pointing;
    m_obs_id     = obs.m_obs_id;
    m_ontime     = obs.m_ontime;
    m_livetime   = obs.m_livetime;
    m_deadc      = obs.m_deadc;
    m_ewidth     = obs.m_ewidth;

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
    const GSkyMap& ref = dynamic_cast<GCOMEventCube*>(m_events)->map();

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
 *     RA_SCZ   - Right Ascension of pointing
 *     DEC_SCZ  - Declination of pointing
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

        // Set pointing information
        double ra_scz  = hdu->real("RA_SCZ");
        double dec_scz = hdu->real("DEC_SCZ");
        m_pointing.radec_deg(ra_scz, dec_scz);

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
