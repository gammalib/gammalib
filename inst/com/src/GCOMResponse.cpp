/***************************************************************************
 *                GCOMResponse.cpp - COMPTEL Response class                *
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
#include "GFitsImageFloat.hpp"
#include "GModelSky.hpp"
#include "GModelSpatialPointSource.hpp"
#include "GModelSpectralConst.hpp"
#include "GCOMResponse.hpp"
#include "GCOMObservation.hpp"
#include "GCOMEventCube.hpp"
#include "GCOMEventBin.hpp"
#include "GCOMInstDir.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_IRF           "GCOMResponse::irf(GEvent&, GPhoton&, GObservation&)"
#define G_IRF_SPATIAL         "GCOMResponse::irf_spatial(GEvent&, GSource&, "\
                                                             "GObservation&)"
#define G_NROI            "GCOMResponse::nroi(GModelSky&, GEnergy&, GTime&, "\
                                                             "GObservation&)"
#define G_EBOUNDS                           "GCOMResponse::ebounds(GEnergy&)"
#define G_IRF_PTSRC     "GCOMResponse::irf_ptsrc(GModelSky&, GObservation&, "\
                                                                  "GMatrix*)"
#define G_IRF_RADIAL   "GCOMResponse::irf_radial(GModelSky&, GObservation&, "\
                                                                  "GMatrix*)"
#define G_IRF_ELLIPTICAL          "GCOMResponse::irf_elliptical(GModelSky&, "\
                                                  " GObservation&, GMatrix*)"
#define G_IRF_DIFFUSE "GCOMResponse::irf_diffuse(GModelSky&, GObservation&, "\
                                                                  "GMatrix*)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */
#define G_USE_DRM_CUBE

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
 * @return Instrument response function (\f$cm^2 sr^{-1}\f$)
 *
 * @exception GException::invalid_argument
 *            Observation is not a COMPTEL observation.
 *            Event is not a COMPTEL event bin.
 * @exception GException::invalid_value
 *            Response not initialised with a valid IAQ
 *
 * Returns the instrument response function for a given observed photon
 * direction as function of the assumed true photon direction. The result
 * is given by
 *
 * \f[
 *    {\tt IRF} = \frac{{\tt IAQ} \times {\tt DRG} \times {\tt DRX}}
 *                     {T \times {\tt TOFCOR}}
 * \f]
 *
 * where
 * - \f${\tt IRF}\f$ is the instrument response function (\f$cm^2 sr^{-1}\f$),
 * - \f${\tt IAQ}\f$ is the COMPTEL response matrix (\f$sr^{-1}\f$),
 * - \f${\tt DRG}\f$ is the geometry factor (probability),
 * - \f${\tt DRX}\f$ is the exposure (\f$cm^2 s\f$),
 * - \f$T\f$ is the ontime (\f$s\f$), and
 * - \f${\tt TOFCOR}\f$ is a correction that accounts for the Time of Flight
 *   selection window.
 *
 * The observed photon direction is spanned by the three values \f$\Chi\f$,
 * \f$\Psi\f$, and \f$\bar{\varphi})\f$. \f$\Chi\f$ and \f$\Psi\f$ is the
 * scatter direction of the event, given in sky coordinates.
 * \f$\bar{\varphi}\f$ is the Compton scatter angle, computed from the
 * energy deposits in the two detector planes.
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

    // Extract COMPTEL event cube
    const GCOMEventCube* cube = dynamic_cast<const GCOMEventCube*>(observation->events());
    if (cube == NULL) {
        std::string cls = std::string(typeid(&cube).name());
        std::string msg = "Event cube of type \""+cls+"\" is  not a COMPTEL "
                          "event cube. Please specify a COMPTEL event cube "
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

    // Throw an exception if COMPTEL response is not set or if
    if (m_iaq.empty()) {
        std::string msg = "COMPTEL response is empty. Please initialise the "
                          "response with an \"IAQ\".";
        throw GException::invalid_value(G_IRF, msg);
    }
    else if (m_phigeo_bin_size == 0.0) {
        std::string msg = "COMPTEL response has a zero Phigeo bin size. "
                          "Please initialise the response with a valid "
                          "\"IAQ\".";
        throw GException::invalid_value(G_IRF, msg);
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

    // Initialise IRF
    double iaq = 0.0;
    double irf = 0.0;

    // Extract IAQ value by linear inter/extrapolation in Phigeo
    if (iphibar < m_phibar_bins) {
        double phirat  = phigeo / m_phigeo_bin_size; // 0.5 at bin centre
        int    iphigeo = int(phirat);                // index into which Phigeo falls
        double eps     = phirat - iphigeo - 0.5;     // 0.0 at bin centre [-0.5, 0.5[
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

    // Continue only if IAQ is positive
    if (iaq > 0.0) {

        // Get DRG value (unit: probability)
        double drg = observation->drg().map()(obsDir.dir(), iphibar);

        // Get DRX value (unit: cm^2 sec)
        double drx = observation->drx().map()(srcDir);

        // Get ontime (unit: s)
        double ontime = observation->ontime();

        // Get ToF correction
        double tofcor = cube->dre().tof_correction();

        // Compute IRF value
        irf = iaq * drg * drx / (ontime * tofcor);

        // Apply deadtime correction
        irf *= obs.deadc(srcTime);

        // Make sure that IRF is positive
        if (irf < 0.0) {
            irf = 0.0;
        }

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

    } // endif: IAQ was positive

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
 * Loads the COMPTEL response with specified name @p rspname.
 *
 * The method first attempts to interpret @p rspname as a file name and to
 * load the corresponding response.
 *
 * If @p rspname is not a FITS file the method searches for an appropriate
 * response in the calibration database. If no appropriate response is found,
 * the method takes the CALDB root path and response name to build the full
 * path to the response file, and tries to load the response from these
 * paths.
 *
 * If also this fails an exception is thrown.
 ***************************************************************************/
void GCOMResponse::load(const std::string& rspname)
{
    // Clear instance but conserve calibration database
    GCaldb caldb = m_caldb;
    clear();
    m_caldb = caldb;

    // Save response name
    m_rspname = rspname;

    // Interpret response name as a FITS file name
    GFilename filename(rspname);

    // If the filename does not exist the try getting the response from the
    // calibration database
    if (!filename.is_fits()) {

        // Get GCaldb response
        filename = m_caldb.filename("","","IAQ","","",rspname);

        // If filename is empty then build filename from CALDB root path
        // and response name
        if (filename.is_empty()) {
            filename = gammalib::filepath(m_caldb.rootdir(), rspname);
            if (!filename.exists()) {
                GFilename testname = filename + ".fits";
                if (testname.exists()) {
                    filename = testname;
                }
            }
        }

    } // endif: response name is not a FITS file

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
 *
 * The IAQ values are divided by the solid angle of a Phigeo bin which is
 * given by
 *
 * \f{eqnarray*}{
 *    \Omega & = & 2 \pi \left[
 *      \left(
 *        1 - \cos \left( \varphi_{\rm geo} +
 *                        \frac{1}{2} \Delta \varphi_{\rm geo} \right)
 *      \right) -
 *      \left(
 *        1 - \cos \left( \varphi_{\rm geo} -
 *                        \frac{1}{2} \Delta \varphi_{\rm geo} \right)
 *      \right) \right] \\
 *     &=& 2 \pi \left[
 *        \cos \left( \varphi_{\rm geo} -
 *                    \frac{1}{2} \Delta \varphi_{\rm geo} \right) -
 *        \cos \left( \varphi_{\rm geo} +
 *                    \frac{1}{2} \Delta \varphi_{\rm geo} \right)
 *      \right] \\
 *     &=& 4 \pi \sin \left( \varphi_{\rm geo} \right)
 *             \sin \left( \frac{1}{2} \Delta \varphi_{\rm geo} \right)
 * \f}
 ***************************************************************************/
void GCOMResponse::read(const GFitsImage& image)
{
    // Continue only if there are two image axis
    if (image.naxis() == 2) {

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
        m_phigeo_min = m_phigeo_ref_value + (1.0-m_phigeo_ref_pixel) *
                       m_phigeo_bin_size;
        m_phibar_min = m_phibar_ref_value + (1.0-m_phibar_ref_pixel) *
                       m_phibar_bin_size;

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

        // Precompute variable for conversion
        double phigeo_min = m_phigeo_min      * gammalib::deg2rad;
        double phigeo_bin = m_phigeo_bin_size * gammalib::deg2rad;
        double omega0     = gammalib::fourpi  * std::sin(0.5 * phigeo_bin);

        // Convert IAQ matrix from probability per Phigeo bin into a
        // probability per steradian
        for (int iphigeo = 0; iphigeo < m_phigeo_bins; ++iphigeo) {
            double phigeo = iphigeo * phigeo_bin + phigeo_min;
            double omega  = omega0  * std::sin(phigeo);
            for (int iphibar = 0; iphibar < m_phibar_bins; ++iphibar) {
                m_iaq[iphigeo+iphibar*m_phigeo_bins] /= omega;
            }
        }

    } // endif: image had two axes

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write COMPTEL response into FITS image.
 *
 * @param[in] image FITS image.
 *
 * Writes the COMPTEL response into an IAQ FITS file.
 ***************************************************************************/
void GCOMResponse::write(GFitsImageFloat& image) const
{
    // Continue only if response is not empty
    if (m_phigeo_bins > 0 && m_phibar_bins > 0) {

        // Initialise image
        image = GFitsImageFloat(m_phigeo_bins, m_phibar_bins);

        // Convert IAQ matrix from probability per steradian into
        //probability per Phigeo bin
        double omega0 = gammalib::fourpi *
                        std::sin(0.5 * m_phigeo_bin_size * gammalib::deg2rad);
        for (int iphigeo = 0; iphigeo < m_phigeo_bins; ++iphigeo) {
            double phigeo = iphigeo * m_phigeo_bin_size + m_phigeo_min;
            double omega  = omega0 * std::sin(phigeo * gammalib::deg2rad);
            for (int iphibar = 0; iphibar < m_phibar_bins; ++iphibar) {
                int i          = iphigeo + iphibar*m_phigeo_bins;
                image(iphigeo, iphibar) = m_iaq[i] * omega;
            }
        }

        // Set header keywords
        image.card("CTYPE1", "Phigeo", "Geometrical scatter angle");
        image.card("CRVAL1", m_phigeo_ref_value,
                   "[deg] Geometrical scatter angle for reference bin");
        image.card("CDELT1", m_phigeo_bin_size,
                   "[deg] Geometrical scatter angle bin size");
        image.card("CRPIX1", m_phigeo_ref_pixel,
                   "Reference bin of geometrical scatter angle");
        image.card("CTYPE2", "Phibar", "Compton scatter angle");
        image.card("CRVAL2", m_phibar_ref_value,
                   "[deg] Compton scatter angle for reference bin");
        image.card("CDELT2", m_phibar_bin_size, "[deg] Compton scatter angle bin size");
        image.card("CRPIX2", m_phibar_ref_pixel,
                   "Reference bin of Compton scatter angle");
        image.card("BUNIT", "Probability per bin", "Unit of bins");
        image.card("TELESCOP", "GRO", "Name of telescope");
        image.card("INSTRUME", "COMPTEL", "Name of instrument");
        image.card("ORIGIN", "GammaLib", "Origin of FITS file");
        image.card("OBSERVER", "Unknown", "Observer that created FITS file");

    } // endif: response was not empty

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


/***********************************************************************//**
 * @brief Return instrument response to point source
 *
 * @param[in] model Sky model.
 * @param[in] obs Observation.
 * @param[out] gradients Gradients matrix.
 * @return Instrument response to point source for all events in
 *         observation (\f$cm^2\f$).
 *
 * @exception GException::invalid_argument
 *            Observation is not a COMPTEL observation.
 *            Event is not a COMPTEL event bin.
 * @exception GException::invalid_value
 *            Response not initialised with a valid IAQ
 *
 * Returns the instrument response to a point source for all events in the
 * observations.
 *
 * @p gradients is an optional matrix where the number of rows corresponds
 * to the number of events in the observation and the number of columns
 * corresponds to the number of spatial model parameters. Since for point
 * sources no gradients are computed, the method does not alter the
 * content of @p gradients.
 ***************************************************************************/
GVector GCOMResponse::irf_ptsrc(const GModelSky&    model,
                                const GObservation& obs,
                                GMatrix*            gradients) const
{
    // Extract COMPTEL observation
    const GCOMObservation* obs_ptr = dynamic_cast<const GCOMObservation*>(&obs);
    if (obs_ptr == NULL) {
        std::string cls = std::string(typeid(&obs).name());
        std::string msg = "Observation of type \""+cls+"\" is not a COMPTEL "
                          "observations. Please specify a COMPTEL observation "
                          "as argument.";
        throw GException::invalid_argument(G_IRF_PTSRC, msg);
    }

    // Extract COMPTEL event cube
    const GCOMEventCube* cube = dynamic_cast<const GCOMEventCube*>(obs_ptr->events());
    if (cube == NULL) {
        std::string msg = "Observation \""+obs.name()+"\" ("+obs.id()+") does "
                          "not contain a COMPTEL event cube. Please specify "
                          "a COMPTEL observation containing and event cube.";
        throw GException::invalid_argument(G_IRF_PTSRC, msg);
    }

    // Throw an exception if COMPTEL response is not set or if
    if (m_iaq.empty()) {
        std::string msg = "COMPTEL response is empty. Please initialise the "
                          "response with an \"IAQ\".";
        throw GException::invalid_value(G_IRF_PTSRC, msg);
    }
    else if (m_phigeo_bin_size == 0.0) {
        std::string msg = "COMPTEL response has a zero Phigeo bin size. "
                          "Please initialise the response with a valid "
                          "\"IAQ\".";
        throw GException::invalid_value(G_IRF_PTSRC, msg);
    }

    // Get number of Chi/Psi pixels, Phibar layers and event bins
    int npix    = cube->naxis(0) * cube->naxis(1);
    int nphibar = cube->naxis(2);
    int nevents = cube->size();

    // Initialise result
    GVector irfs(nevents);

    // Get point source direction
    GSkyDir srcDir =
    static_cast<const GModelSpatialPointSource*>(model.spatial())->dir();

    // Get IAQ normalisation (cm2): DRX (cm2 s) * DEADC / ONTIME (s)
    double iaq_norm = obs_ptr->drx().map()(srcDir) * obs_ptr->deadc() /
                      (obs_ptr->ontime() * cube->dre().tof_correction());

    // Get pointer to DRG pixels
    const double* drg = obs_ptr->drg().map().pixels();

    // Loop over Chi and Psi
    for (int ipix = 0; ipix < npix; ++ipix) {

        // Get pointer to event bin
        const GCOMEventBin* bin = (*cube)[ipix];

        // Get reference to instrument direction
        const GCOMInstDir& obsDir = bin->dir();

        // Get reference to Phigeo sky direction
        const GSkyDir& phigeoDir = obsDir.dir();

        // Compute angle between true photon arrival direction and scatter
        // direction (Chi,Psi)
        double phigeo = srcDir.dist_deg(phigeoDir);

        // Compute interpolation factors
        double phirat  = phigeo / m_phigeo_bin_size; // 0.5 at bin centre
        int    iphigeo = int(phirat);                // index into which Phigeo falls
        double eps     = phirat - iphigeo - 0.5;     // 0.0 at bin centre [-0.5, 0.5[

        // If Phigeo is in range then compute the IRF value for all
        // Phibar layers
        if (iphigeo < m_phigeo_bins) {

            // Loop over Phibar
            for (int iphibar = 0; iphibar < nphibar; ++iphibar) {

                // Get IAQ index
                int i = iphibar * m_phigeo_bins + iphigeo;

                // Initialise IAQ
                double iaq = 0.0;

                // Compute IAQ
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

                // Continue only if IAQ is positive
                if (iaq > 0.0) {

                    // Get DRI index
                    int idri = ipix + iphibar * npix;

                    // Compute IRF value
                    double irf = iaq * drg[idri] * iaq_norm;

                    // Make sure that IRF is positive
                    if (irf < 0.0) {
                        irf = 0.0;
                    }

                    // Store IRF value
                    irfs[idri] = irf;

                } // endif: IAQ was positive

            } // endfor: looped over Phibar

        } // endif: Phigeo was in valid range

    } // endfor: looped over Chi and Psi

    // Return IRF vector
    return irfs;
}


/***********************************************************************//**
 * @brief Return instrument response to radial source
 *
 * @param[in] model Sky model.
 * @param[in] obs Observation.
 * @param[out] gradients Gradients matrix.
 * @return Instrument response to radial source  for all events in
 *         observation (\f$cm^2\f$).
 *
 * @todo Implement method.
 ***************************************************************************/
GVector GCOMResponse::irf_radial(const GModelSky&    model,
                                 const GObservation& obs,
                                 GMatrix*            gradients) const
{
    // Throw exception
    std::string msg = "Response computation not yet implemented "
                      "for spatial model type \""+
                      model.spatial()->type()+"\".";
    throw GException::feature_not_implemented(G_IRF_RADIAL, msg);

    // Get number of events
    int nevents = obs.events()->size();

    // Initialise result
    GVector irfs(nevents);

    // Return IRF value
    return irfs;
}


/***********************************************************************//**
 * @brief Return instrument response to elliptical source
 *
 * @param[in] model Sky model.
 * @param[in] obs Observation.
 * @param[out] gradients Gradients matrix.
 * @return Instrument response to elliptical source for all events in
 *         observation (\f$cm^2\f$).
 *
 * @todo Implement method.
 ***************************************************************************/
GVector GCOMResponse::irf_elliptical(const GModelSky&    model,
                                     const GObservation& obs,
                                     GMatrix*            gradients) const
{
    // Throw exception
    std::string msg = "Response computation not yet implemented "
                      "for spatial model type \""+
                      model.spatial()->type()+"\".";
    throw GException::feature_not_implemented(G_IRF_ELLIPTICAL, msg);

    // Get number of events
    int nevents = obs.events()->size();

    // Initialise result
    GVector irfs(nevents);

    // Return IRF value
    return irfs;
}


/***********************************************************************//**
 * @brief Return instrument response to diffuse source
 *
 * @param[in] model Sky model.
 * @param[in] obs Observation.
 * @param[out] gradients Gradients matrix.
 * @return Instrument response to diffuse source for all events in
 *         observation (\f$cm^2\f$).
 *
 * @exception GException::invalid_argument
 *            Observation is not a COMPTEL observation.
 * @exception GException::invalid_value
 *            Response not initialised with a valid IAQ
 *
 * Returns the instrument response to a diffuse source for all events in
 * the observations. The diffuse source may be energy dependent.
 *
 * The computation is done by integrating the diffuse model for each pixel
 * in Chi and Psi over a circular region centred on the Chi/Psi pixel with
 * a radius equal to the maximum Phigeo value. The radial integration is
 * done by looping over all Phigeo bins of the response. For each Phigeo
 * value, the azimuthal integration is done by stepping with an angular step
 * size that corresponds to the Phigeo step size (which typically is 1 deg).
 *
 * @p gradients is an optional matrix where the number of rows corresponds
 * to the number of events in the observation and the number of columns
 * corresponds to the number of spatial model parameters. Since for point
 * sources no gradients are computed, the method does not alter the
 * content of @p gradients.
 ***************************************************************************/
GVector GCOMResponse::irf_diffuse(const GModelSky&    model,
                                  const GObservation& obs,
                                  GMatrix*            gradients) const
{
    // Extract COMPTEL observation
    const GCOMObservation* obs_ptr = dynamic_cast<const GCOMObservation*>(&obs);
    if (obs_ptr == NULL) {
        std::string cls = std::string(typeid(&obs).name());
        std::string msg = "Observation of type \""+cls+"\" is not a COMPTEL "
                          "observations. Please specify a COMPTEL observation "
                          "as argument.";
        throw GException::invalid_argument(G_IRF_DIFFUSE, msg);
    }

    // Extract COMPTEL event cube
    const GCOMEventCube* cube = dynamic_cast<const GCOMEventCube*>(obs_ptr->events());
    if (cube == NULL) {
        std::string msg = "Observation \""+obs.name()+"\" ("+obs.id()+") does "
                          "not contain a COMPTEL event cube. Please specify "
                          "a COMPTEL observation containing and event cube.";
        throw GException::invalid_argument(G_IRF_DIFFUSE, msg);
    }

    // Throw an exception if COMPTEL response is not set or if
    if (m_iaq.empty()) {
        std::string msg = "COMPTEL response is empty. Please initialise the "
                          "response with an \"IAQ\".";
        throw GException::invalid_value(G_IRF_DIFFUSE, msg);
    }
    else if (m_phigeo_bin_size == 0.0) {
        std::string msg = "COMPTEL response has a zero Phigeo bin size. "
                          "Please initialise the response with a valid "
                          "\"IAQ\".";
        throw GException::invalid_value(G_IRF_DIFFUSE, msg);
    }

    // Get number of Chi/Psi pixels, Phibar layers and event bins
    int npix    = cube->naxis(0) * cube->naxis(1);
    int nphibar = cube->naxis(2);
    int nevents = cube->size();

    // Initialise result
    GVector irfs(nevents);

    // Initialise some variables
    double         phigeo_min = m_phigeo_min      * gammalib::deg2rad;
    double         phigeo_bin = m_phigeo_bin_size * gammalib::deg2rad;
    double         omega0     = 2.0 * std::sin(0.5 * phigeo_bin);
    const GSkyMap& drx        = obs_ptr->drx().map();
    const double*  drg        = obs_ptr->drg().map().pixels();

    // Compute IAQ normalisation (1/s): DEADC / ONTIME (s)
    double iaq_norm = obs_ptr->deadc() /
                      (obs_ptr->ontime() * cube->dre().tof_correction());

    // Loop over Chi and Psi
    for (int ipix = 0; ipix < npix; ++ipix) {

        // Get pointer to event bin
        const GCOMEventBin* bin = (*cube)[ipix];

        // Get reference to instrument direction
        const GCOMInstDir& obsDir = bin->dir();

        // Get reference to Phigeo sky direction
        const GSkyDir& phigeoDir = obsDir.dir();

        // Loop over Phigeo
        for (int iphigeo = 0; iphigeo < m_phigeo_bins; ++iphigeo) {

            // Determine Phigeo value in radians
            double phigeo     = phigeo_min + iphigeo * phigeo_bin;
            double sin_phigeo = std::sin(phigeo);

            // Determine number of azimuthal integration steps and step size
            // in radians
            double length = gammalib::twopi * sin_phigeo;
            int    naz    = int(length / phigeo_bin + 0.5);
            if (naz < 2) {
                naz = 2;
            }
            double az_step = gammalib::twopi / double(naz);

            // Computes solid angle of integration bin multiplied by Jaccobian
            double omega = omega0 * az_step * sin_phigeo;

            // Loop over azimuth angle
            double az = 0.0;
            for (int iaz = 0; iaz < naz; ++iaz, az += az_step) {

                // Get sky direction
                GSkyDir skyDir = phigeoDir;
                skyDir.rotate(az, phigeo);

                // Fall through if sky direction is not contained in DRX
                if (!drx.contains(skyDir)) {
                    continue;
                }

                // Set photon
                GPhoton photon(skyDir, bin->energy(), bin->time());

                // Get model sky intensity for photon (unit: sr^-1)
                double intensity = model.spatial()->eval(photon);

                // Fall through if intensity is zero
                if (intensity == 0.0) {
                    continue;
                }

                // Multiply intensity by DRX value (unit: cm^2 s sr^-1)
                intensity *= drx(skyDir);

                // Multiply intensity by solid angle (unit: cm^2 s)
                intensity *= omega;

                // Loop over Phibar
                for (int iphibar = 0; iphibar < nphibar; ++iphibar) {

                    // Get IAQ index
                    int i = iphibar * m_phigeo_bins + iphigeo;

                    // Get IAQ value
                    double iaq = m_iaq[i];

                    // Fall through if IAQ is not positive
                    if (iaq <= 0.0) {
                        continue;
                    }

                    // Get DRI index
                    int idri = ipix + iphibar * npix;

                    // Compute IRF value (unit: cm^2)
                    double irf = iaq * drg[idri] * iaq_norm * intensity;

                    // Add IRF value if it is positive
                    if (irf > 0.0) {
                        irfs[idri] += irf;
                    }

                } // endfor: looped over Phibar

            } // endfor: looped over azimuth angle around locus of IAQ

        } // endfor: looped over Phigeo angles

    } // endfor: looped over Chi and Psi pixels

    // Return IRF vector
    return irfs;
}
