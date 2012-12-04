/***************************************************************************
 *               GCOMResponse.cpp  -  COMPTEL Response class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012 by Juergen Knoedlseder                              *
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
#include "GTools.hpp"
#include "GFits.hpp"
#include "GCaldb.hpp"
#include "GCOMResponse.hpp"
#include "GCOMObservation.hpp"
#include "GCOMInstDir.hpp"
#include "GCOMException.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_IRF         "GCOMResponse::irf(GInstDir&,GEnergy&,GTime&,GSkyDir&,"\
                                             "GEnergy&,GTime&,GObservation&)"
#define G_NPRED               "GCOMResponse::npred(GSkyDir&,GEnergy&,GTime&,"\
                                                             "GObservation&)"

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
 * @param[in] iaqname IAQ file name.
 * @param[in] caldb Calibration database path (defaults to "").
 *
 * Create COMPTEL response by loading an IAQ file.
 *
 * If an empty string is passed as calibration database path, the method will
 * use the CALDB environment variable to determine the calibration database
 * path. This is done in the method GCOMResponse::caldb which makes use of a
 * GCaldb object to locate the calibration database.
 ***************************************************************************/
GCOMResponse::GCOMResponse(const std::string& iaqname,
                           const std::string& caldb) : GResponse()
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
 * Creates a clone (deep copy) of a COMPTEL response object.
 ***************************************************************************/
GCOMResponse* GCOMResponse::clone(void) const
{
    return new GCOMResponse(*this);
}


/***********************************************************************//**
 * @brief Return value of instrument response function
 *
 * @param[in] obsDir Observed photon direction.
 * @param[in] obsEng Observed energy of photon (not used).
 * @param[in] obsTime Observed photon arrival time (not used).
 * @param[in] srcDir True photon arrival direction.
 * @param[in] srcEng True energy of photon (not used).
 * @param[in] srcTime True photon arrival time (not used).
 * @param[in] obs Observation.
 * @return Instrument response function (cm2 sr-1)
 *
 * @exception GCOMException::bad_observation_type
 *            Observation is not a COMPTEL observation.
 * @exception GCOMException::bad_instdir_type
 *            Instrument direction is not a COMPTEL instrument direction.
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
double GCOMResponse::irf(const GInstDir&     obsDir,
                         const GEnergy&      obsEng,
                         const GTime&        obsTime,
                         const GSkyDir&      srcDir,
                         const GEnergy&      srcEng,
                         const GTime&        srcTime,
                         const GObservation& obs) const
{
    // Extract COMPTEL observation
    const GCOMObservation* observation = dynamic_cast<const GCOMObservation*>(&obs);
    if (observation == NULL) {
        throw GCOMException::bad_observation_type(G_IRF);
    }

    // Extract COMPTEL instrument direction
    const GCOMInstDir* dir = dynamic_cast<const GCOMInstDir*>(&obsDir);
    if (dir == NULL) {
        throw GCOMException::bad_instdir_type(G_IRF);
    }

    // Compute angle between true photon arrival direction and scatter
    // direction (Chi,Psi)
    double phigeo = srcDir.dist_deg(dir->skydir());

    // Compute scatter angle index
    int iphibar = int(dir->phi() / m_phibar_bin_size);

    // Extract IAQ value by linear inter/extrapolation in Phigeo
    double phirat  = phigeo / m_phigeo_bin_size; // 0.5 at bin centre
    int    iphigeo = int(phirat);                // index into which Phigeo falls
    double eps     = phirat - iphigeo - 0.5;     // 0.0 at bin centre
    double iaq     = 0.0;
    if (iphigeo < m_phigeo_bins) {
        int i = iphibar * m_phigeo_bins + iphigeo;
        if (eps < 0.0 && iphigeo > 0) { // interpolate towards left
            iaq = (1.0 + eps) * m_iaq[i] - eps * m_iaq[i-1];
        }
        else {                          // interpolate towards right
            iaq = (1.0 - eps) * m_iaq[i] + eps * m_iaq[i+1];
        }
    }

    // Get DRG value (units: cm2)
    double drg = observation->drg()(dir->skydir(), iphibar);

    // Get DRX value (units: sec)
    double drx = observation->drx()(srcDir);

    // Get ontime
    double ontime = observation->ontime(); // sec

    // Compute IRF value
    double irf = iaq * drg * drx / ontime;

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (isnotanumber(irf) || isinfinite(irf)) {
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
 * @brief Return spatial integral of point spread function
 *
 * @param[in] srcDir True photon arrival direction.
 * @param[in] srcEng True energy of photon.
 * @param[in] srcTime True photon arrival time.
 * @param[in] obs Observation.
 *
 * @todo Implement method (is maybe not really needed)
 ***************************************************************************/
double GCOMResponse::npred(const GSkyDir&      srcDir,
                           const GEnergy&      srcEng,
                           const GTime&        srcTime,
                           const GObservation& obs) const
{
    // Set dummp Npred value
    double npred = 1.0;

    // Compile option: Check for NaN/Inf
    #if defined(G_NAN_CHECK)
    if (isnotanumber(npred) || isinfinite(npred)) {
        std::cout << "*** ERROR: GCOMResponse::npred:";
        std::cout << " NaN/Inf encountered";
        std::cout << " (";
        std::cout << "npred=" << npred;
        std::cout << ")";
        std::cout << std::endl;
    }
    #endif

    // Return Npred
    return npred;
}


/***********************************************************************//**
 * @brief Set path to the calibration database
 *
 * @param[in] caldb Path to calibration database
 *
 * This method stores the CALDB root directory as the path to the COMPTEL
 * calibration database. Once the final calibration format has been decided
 * on, this method should be adjusted.
 ***************************************************************************/
void GCOMResponse::caldb(const std::string& caldb)
{
    // Allocate calibration database
    GCaldb db(caldb);

    // Store the path to the calibration database
    m_caldb = db.dir();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return path to the calibration database
 ***************************************************************************/
std::string GCOMResponse::caldb(void) const
{
    // Return
    return m_caldb;
}


/***********************************************************************//**
 * @brief Load COMPTEL response.
 *
 * @param[in] iaqname Name of IAQ file.
 *
 * Load COMPTEL response from IAQ file.
 ***************************************************************************/
void GCOMResponse::load(const std::string& iaqname)
{
    // Save calibration database name
    std::string caldb = m_caldb;

    // Clear instance
    clear();

    // Restore calibration database name
    m_caldb = caldb;

    // Save IAQ name
    m_iaqname = iaqname;

    // Build filename
    std::string filename = m_caldb + "/" + m_iaqname;

    // Open FITS file
    GFits file(filename);

    // Get IAQ image
    GFitsImage* iaq = file.image(0);

    // Read IAQ
    read_iaq(iaq);

    // Close ARF FITS file
    file.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read COMPTEL response from FITS image.
 *
 * @param[in] hdu FITS Image HDU.
 *
 * Read the COMPTEL response from IAQ FITS file and convert the IAQ values
 * into a probability per steradian.
 ***************************************************************************/
void GCOMResponse::read_iaq(const GFitsImage* hdu)
{
    // Continue only if header is valid
    if (hdu != NULL) {

        // Store IAQ dimensions
        m_phigeo_bins = hdu->naxes(0);
        m_phibar_bins = hdu->naxes(1);

        // Store IAQ axes definitions
        m_phigeo_ref_value = hdu->real("CRVAL1");
        m_phigeo_ref_pixel = hdu->real("CRPIX1");
        m_phigeo_bin_size  = hdu->real("CDELT1");
        m_phibar_ref_value = hdu->real("CRVAL2");
        m_phibar_ref_pixel = hdu->real("CRPIX2");
        m_phibar_bin_size  = hdu->real("CDELT2");

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
                m_iaq[i] = hdu->pixel(i);
            }

        } // endif: size was positive

        // Convert IAQ matrix from probability per Phigeo bin into a
        // probability per steradian
        double omega0 = fourpi * std::sin(0.5 * m_phigeo_bin_size * deg2rad);
        for (int iphigeo = 0; iphigeo < m_phigeo_bins; ++iphigeo) {
            double phigeo = iphigeo * m_phigeo_bin_size + m_phigeo_min;
            double omega  = omega0 * std::sin(phigeo * deg2rad);
            for (int iphibar = 0; iphibar < m_phibar_bins; ++iphibar) {
                m_iaq[iphigeo+iphibar*m_phigeo_bins] /= omega;
            }
        }

    } // endif: HDU was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print COMPTEL response information
 ***************************************************************************/
std::string GCOMResponse::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GCOMResponse ===");
    result.append("\n"+parformat("Calibration database")+m_caldb);
    result.append("\n"+parformat("IAQ file name")+m_iaqname);
    result.append("\n"+parformat("Number of Phigeo bins")+str(m_phigeo_bins));
    result.append("\n"+parformat("Number of Phibar bins")+str(m_phibar_bins));
    result.append("\n"+parformat("Phigeo reference value")+str(m_phigeo_ref_value)+" deg");
    result.append("\n"+parformat("Phigeo reference pixel")+str(m_phigeo_ref_pixel));
    result.append("\n"+parformat("Phigeo bin size")+str(m_phigeo_bin_size)+" deg");
    result.append("\n"+parformat("Phigeo first bin value")+str(m_phigeo_min)+" deg");
    result.append("\n"+parformat("Phibar reference value")+str(m_phibar_ref_value)+" deg");
    result.append("\n"+parformat("Phibar reference pixel")+str(m_phibar_ref_pixel));
    result.append("\n"+parformat("Phibar bin size")+str(m_phibar_bin_size)+" deg");
    result.append("\n"+parformat("Phibar first bin value")+str(m_phibar_min)+" deg");

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
    m_iaqname.clear();
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
    m_iaqname          = rsp.m_iaqname;
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
