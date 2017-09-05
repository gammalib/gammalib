/***************************************************************************
 *                GCOMResponse.cpp - COMPTEL Response class                *
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
 * @file GCOMResponse.cpp
 * @brief COMPTEL response class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <string>
#include <typeinfo>
#include "GTools.hpp"
#include "GMath.hpp"
#include "GFits.hpp"
#include "GCaldb.hpp"
#include "GEvent.hpp"
#include "GPhoton.hpp"
#include "GSource.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GObservation.hpp"
#include "GFitsImage.hpp"
#include "GModelSky.hpp"
#include "GModelSpatialPointSource.hpp"
#include "GCOMResponse.hpp"
#include "GCOMObservation.hpp"
#include "GCOMEventBin.hpp"
#include "GCOMInstDir.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_IRF           "GCOMResponse::irf(GEvent&, GSource&, GObservation&)"
#define G_NROI            "GCOMResponse::nroi(GModelSky&, GEnergy&, GTime&, "\
                                                             "GObservation&)"
#define G_EBOUNDS                           "GCOMResponse::ebounds(GEnergy&)"

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
 *
 * Creates an empty COMPTEL response.
 ***************************************************************************/
GCOMResponse::GCOMResponse(void) : GResponse()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] rsp COM response.
 **************************************************************************/
GCOMResponse::GCOMResponse(const GCOMResponse& rsp) : GResponse(rsp)
{
    // Initialise members
    init_members();

    // Copy members
    copy_members(rsp);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Response constructor
 *
 * @param[in] caldb Calibration database.
 * @param[in] iaqname IAQ file name.
 *
 * Create COMPTEL response by loading an IAQ file from a calibration
 * database.
 ***************************************************************************/
GCOMResponse::GCOMResponse(const GCaldb&      caldb,
                           const std::string& iaqname) : GResponse()
{
    // Initialise members
    init_members();

    // Set calibration database
    this->caldb(caldb);

    // Load IRF
    this->load(iaqname);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 *
 * Destroys instance of COMPTEL response object.
 ***************************************************************************/
GCOMResponse::~GCOMResponse(void)
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
 * @param[in] rsp COMPTEL response.
 * @return COMPTEL response.
 *
 * Assigns COMPTEL response object to another COMPTEL response object. The
 * assignment performs a deep copy of all information, hence the original
 * object from which the assignment has been performed can be destroyed after
 * this operation without any loss of information.
 ***************************************************************************/
GCOMResponse& GCOMResponse::operator=(const GCOMResponse& rsp)
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
 * @brief Clear instance
 *
 * Clears COMPTEL response object by resetting all members to an initial
 * state. Any information that was present in the object before will be lost.
 ***************************************************************************/
void GCOMResponse::clear(void)
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
 *
 * @return Pointer to deep copy of COMPTEL response.
 ***************************************************************************/
GCOMResponse* GCOMResponse::clone(void) const
{
    return new GCOMResponse(*this);
}


/***********************************************************************//**
 * @brief Return value of instrument response function
 *
 * @param[in] event Observed event.
 * @param[in] photon Incident photon.
 * @param[in] obs Observation.
 * @return Instrument response function (cm2 sr-1)
 *
 * @exception GException::invalid_argument
 *            Observation is not a COMPTEL observation.
 *            Event is not a COMPTEL event bin.
 *
 * Returns the instrument response function for a given observed photon
 * direction as function of the assumed true photon direction. The result
 * is given by
 * \f[IRF = \frac{IAQ \times DRG \times DRX}{ontime \times ewidth}\f]
 * where
 * \f$IRF\f$ is the instrument response function,
 * \f$IAQ\f$ is the COMPTEL response matrix (sr-1),
 * \f$DRG\f$ is the geometry factor (cm2),
 * \f$DRX\f$ is the exposure (s),
 * \f$ontime\f$ is the ontime (s), and
 * \f$ewidth\f$ is the energy width (MeV).
 *
 * The observed photon direction is spanned by the 3 values (Chi,Psi,Phibar).
 * (Chi,Psi) is the scatter direction of the event, given in sky coordinates.
 * Phibar is the Compton scatter angle, computed from the energy deposits.
 ***************************************************************************/
double GCOMResponse::irf(const GEvent&       event,
                         const GPhoton&      photon,
                         const GObservation& obs) const
{
    // Extract COMPTEL observation
    const GCOMObservation* observation = dynamic_cast<const GCOMObservation*>(&obs);
    if (observation == NULL) {
        std::string cls = std::string(typeid(&obs).name());
        std::string msg = "Observation of type \""+cls+"\" is not a COMPTEL "
                          "observations. Please specify a COMPTEL observation "
                          "as argument.";
        throw GException::invalid_argument(G_IRF, msg);
    }

    // Extract COMPTEL event bin
    const GCOMEventBin* bin = dynamic_cast<const GCOMEventBin*>(&event);
    if (bin == NULL) {
        std::string cls = std::string(typeid(&event).name());
        std::string msg = "Event of type \""+cls+"\" is  not a COMPTEL event. "
                          "Please specify a COMPTEL event as argument.";
        throw GException::invalid_argument(G_IRF, msg);
    }

    // Extract event parameters
    const GCOMInstDir& obsDir = bin->dir();

    // Extract photon parameters
    const GSkyDir& srcDir  = photon.dir();
    const GTime&   srcTime = photon.time();

    // Compute angle between true photon arrival direction and scatter
    // direction (Chi,Psi)
    double phigeo = srcDir.dist_deg(obsDir.dir());

    // Compute scatter angle index
    int iphibar = int(obsDir.phibar() / m_phibar_bin_size);

    // Extract IAQ value by linear inter/extrapolation in Phigeo
    double iaq = 0.0;
    if (iphibar < m_phibar_bins) {
        double phirat  = phigeo / m_phigeo_bin_size; // 0.5 at bin centre
        int    iphigeo = int(phirat);                // index into which Phigeo falls
        double eps     = phirat - iphigeo - 0.5;     // 0.0 at bin centre
        if (iphigeo < m_phigeo_bins) {
            int i = iphibar * m_phigeo_bins + iphigeo;
            if (eps < 0.0) { // interpolate towards left
                if (iphigeo > 0) {
                    iaq = (1.0 + eps) * m_iaq[i] - eps * m_iaq[i-1];
                }
                else {
                    iaq = (1.0 - eps) * m_iaq[i] + eps * m_iaq[i+1];
                }
            }
            else {           // interpolate towards right
                if (iphigeo < m_phigeo_bins-1) {
                    iaq = (1.0 - eps) * m_iaq[i] + eps * m_iaq[i+1];
                }
                else {
                    iaq = (1.0 + eps) * m_iaq[i] - eps * m_iaq[i-1];
                }
            }
        }
    }

    // Get DRG value (units: cm2)
    double drg = observation->drg()(obsDir.dir(), iphibar);

    // Get DRX value (units: sec)
    double drx = observation->drx()(srcDir);

    // Get ontime
    double ontime = observation->ontime(); // sec

    // Compute IRF value
    double irf = iaq * drg * drx / ontime;

    // Apply deadtime correction
    irf *= obs.deadc(srcTime);

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (gammalib::is_notanumber(irf) || gammalib::is_infinite(irf)) {
        std::cout << "*** ERROR: GCOMResponse::irf:";
        std::cout << " NaN/Inf encountered";
        std::cout << " (irf=" << irf;
        std::cout << ", iaq=" << iaq;
        std::cout << ", drg=" << drg;
        std::cout << ", drx=" << drx;
        std::cout << ")";
        std::cout << std::endl;
    }
    #endif

    // Return IRF value
    return irf;
}


/***********************************************************************//**
 * @brief Return instrument response
 *
 * @param[in] event Event.
 * @param[in] source Source.
 * @param[in] obs Observation.
 * @return Instrument response.
 *
 * Returns the instrument response for a given event, source and observation.
 ***************************************************************************/
double GCOMResponse::irf(const GEvent&       event,
                         const GSource&      source,
                         const GObservation& obs) const
{
    // Initialise IRF value
    double irf = 0.0;

    // Select IRF depending on the spatial model type
    switch (source.model()->code()) {
        case GMODEL_SPATIAL_POINT_SOURCE:
            {
            const GModelSpatialPointSource* src =
                  static_cast<const GModelSpatialPointSource*>(source.model());
            GPhoton photon(src->dir(), source.energy(), source.time());
            irf = this->irf(event, photon, obs);
            }
            break;
        case GMODEL_SPATIAL_RADIAL:
        case GMODEL_SPATIAL_ELLIPTICAL:
        case GMODEL_SPATIAL_DIFFUSE:
            {
            std::string msg = "Response computation not yet implemented for "
                              "spatial model type \""+source.model()->type()+"\".";
            throw GException::feature_not_implemented(G_IRF, msg);
            }
            break;
        default:
            break;
    }

    // Return IRF value
    return irf;
}


/***********************************************************************//**
 * @brief Return integral of event probability for a given sky model over ROI
 *
 * @param[in] model Sky model.
 * @param[in] obsEng Observed photon energy.
 * @param[in] obsTime Observed photon arrival time.
 * @param[in] obs Observation.
 * @return 0.0
 *
 * @exception GException::feature_not_implemented
 *            Method is not implemented.
 ***************************************************************************/
double GCOMResponse::nroi(const GModelSky&    model,
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
GEbounds GCOMResponse::ebounds(const GEnergy& obsEnergy) const
{
    // Initialise an empty boundary object
    GEbounds ebounds;

    // Throw an exception
    std::string msg = "Energy dispersion not implemented.";
    throw GException::feature_not_implemented(G_EBOUNDS, msg);

    // Return energy boundaries
    return ebounds;
}


/***********************************************************************//**
 * @brief Load COMPTEL response.
 *
 * @param[in] rspname COMPTEL response name.
 *
 * Loads the COMPTEL response with specified name @p rspname. The method first
 * searchs for an appropriate response in the calibration database. If no
 * appropriate response is found, the method takes the database root path
 * and response name to build the full path to the response file, and tries
 * to load the response from these paths.
 ***************************************************************************/
void GCOMResponse::load(const std::string& rspname)
{
    // Clear instance but conserve calibration database
    GCaldb caldb = m_caldb;
    clear();
    m_caldb = caldb;

    // Save response name
    m_rspname = rspname;

    // First attempt reading the response using the GCaldb interface
    GFilename filename = m_caldb.filename("","","IAQ","","",rspname);

    // If filename is empty then build filename from CALDB root path and
    // response name
    if (filename.is_empty()) {
        filename = gammalib::filepath(m_caldb.rootdir(), rspname);
        if (!filename.exists()) {
            GFilename testname = filename + ".fits";
            if (testname.exists()) {
                filename = testname;
            }
        }
    }

    // Open FITS file
    GFits fits(filename);

    // Get IAQ image
    const GFitsImage& iaq = *fits.image(0);

    // Read IAQ
    read(iaq);

    // Close IAQ FITS file
    fits.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read COMPTEL response from FITS image.
 *
 * @param[in] image FITS image.
 *
 * Read the COMPTEL response from IAQ FITS file and convert the IAQ values
 * into a probability per steradian.
 ***************************************************************************/
void GCOMResponse::read(const GFitsImage& image)
{
    // Store IAQ dimensions
    m_phigeo_bins = image.naxes(0);
    m_phibar_bins = image.naxes(1);

    // Store IAQ axes definitions
    m_phigeo_ref_value = image.real("CRVAL1");
    m_phigeo_ref_pixel = image.real("CRPIX1");
    m_phigeo_bin_size  = image.real("CDELT1");
    m_phibar_ref_value = image.real("CRVAL2");
    m_phibar_ref_pixel = image.real("CRPIX2");
    m_phibar_bin_size  = image.real("CDELT2");

    // Get axes minima (values of first bin)
    m_phigeo_min = m_phigeo_ref_value + (1.0-m_phigeo_ref_pixel) * m_phigeo_bin_size;
    m_phibar_min = m_phibar_ref_value + (1.0-m_phibar_ref_pixel) * m_phibar_bin_size;

    // Compute IAQ size. Continue only if size is positive
    int size = m_phigeo_bins * m_phibar_bins;
    if (size > 0) {

        // Allocate memory for IAQ
        m_iaq.assign(size, 0.0);

        // Copy over IAQ values
        for (int i = 0; i < size; ++i) {
            m_iaq[i] = image.pixel(i);
        }

    } // endif: size was positive

    // Convert IAQ matrix from probability per Phigeo bin into a
    // probability per steradian
    double omega0 = gammalib::fourpi *
                    std::sin(0.5 * m_phigeo_bin_size * gammalib::deg2rad);
    for (int iphigeo = 0; iphigeo < m_phigeo_bins; ++iphigeo) {
        double phigeo = iphigeo * m_phigeo_bin_size + m_phigeo_min;
        double omega  = omega0 * std::sin(phigeo * gammalib::deg2rad);
        for (int iphibar = 0; iphibar < m_phibar_bins; ++iphibar) {
            m_iaq[iphigeo+iphibar*m_phigeo_bins] /= omega;
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print COMPTEL response information
 *
 * @param[in] chatter Chattiness.
 * @return String containing COMPTEL response information.
 ***************************************************************************/
std::string GCOMResponse::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCOMResponse ===");

        // Append response name
        result.append("\n"+gammalib::parformat("Response name")+m_rspname);

        // EXPLICIT: Append detailed information
        if (chatter >= EXPLICIT) {

            // Append information
            result.append("\n"+gammalib::parformat("Number of Phigeo bins"));
            result.append(gammalib::str(m_phigeo_bins));
            result.append("\n"+gammalib::parformat("Number of Phibar bins"));
            result.append(gammalib::str(m_phibar_bins));
            result.append("\n"+gammalib::parformat("Phigeo reference value"));
            result.append(gammalib::str(m_phigeo_ref_value)+" deg");
            result.append("\n"+gammalib::parformat("Phigeo reference pixel"));
            result.append(gammalib::str(m_phigeo_ref_pixel));
            result.append("\n"+gammalib::parformat("Phigeo bin size"));
            result.append(gammalib::str(m_phigeo_bin_size)+" deg");
            result.append("\n"+gammalib::parformat("Phigeo first bin value"));
            result.append(gammalib::str(m_phigeo_min)+" deg");
            result.append("\n"+gammalib::parformat("Phibar reference value"));
            result.append(gammalib::str(m_phibar_ref_value)+" deg");
            result.append("\n"+gammalib::parformat("Phibar reference pixel"));
            result.append(gammalib::str(m_phibar_ref_pixel));
            result.append("\n"+gammalib::parformat("Phibar bin size"));
            result.append(gammalib::str(m_phibar_bin_size)+" deg");
            result.append("\n"+gammalib::parformat("Phibar first bin value"));
            result.append(gammalib::str(m_phibar_min)+" deg");

        }

        // VERBOSE: Append calibration database
        if (chatter == VERBOSE) {
            result.append("\n"+m_caldb.print(chatter));
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
void GCOMResponse::init_members(void)
{
    // Initialise members
    m_caldb.clear();
    m_rspname.clear();
    m_iaq.clear();
    m_phigeo_bins      = 0;
    m_phibar_bins      = 0;
    m_phigeo_ref_value = 0.0;
    m_phigeo_ref_pixel = 0.0;
    m_phigeo_bin_size  = 0.0;
    m_phigeo_min       = 0.0;
    m_phibar_ref_value = 0.0;
    m_phibar_ref_pixel = 0.0;
    m_phibar_bin_size  = 0.0;
    m_phibar_min       = 0.0;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] rsp COMPTEL response.
 ***************************************************************************/
void GCOMResponse::copy_members(const GCOMResponse& rsp)
{
    // Copy attributes
    m_caldb            = rsp.m_caldb;
    m_rspname          = rsp.m_rspname;
    m_iaq              = rsp.m_iaq;
    m_phigeo_bins      = rsp.m_phigeo_bins;
    m_phibar_bins      = rsp.m_phibar_bins;
    m_phigeo_ref_value = rsp.m_phigeo_ref_value;
    m_phigeo_ref_pixel = rsp.m_phigeo_ref_pixel;
    m_phigeo_bin_size  = rsp.m_phigeo_bin_size;
    m_phigeo_min       = rsp.m_phigeo_min;
    m_phibar_ref_value = rsp.m_phibar_ref_value;
    m_phibar_ref_pixel = rsp.m_phibar_ref_pixel;
    m_phibar_bin_size  = rsp.m_phibar_bin_size;
    m_phibar_min       = rsp.m_phibar_min;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCOMResponse::free_members(void)
{
    // Return
    return;
}
