/***************************************************************************
 *              GSPIResponse.cpp - INTEGRAL/SPI response class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2020 by Juergen Knoedlseder                              *
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
 * @file GSPIResponse.hpp
 * @brief INTEGRAL/SPI instrument response function class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <string>
#include <typeinfo>
#include <algorithm>
#include "GException.hpp"
#include "GTools.hpp"
#include "GEbounds.hpp"
#include "GEnergies.hpp"
#include "GNodeArray.hpp"
#include "GFits.hpp"
#include "GFitsBinTable.hpp"
#include "GFitsTableShortCol.hpp"
#include "GFitsTableDoubleCol.hpp"
#include "GSource.hpp"
#include "GModelSpatialPointSource.hpp"
#include "GSPITools.hpp"
#include "GSPIResponse.hpp"
#include "GSPIObservation.hpp"
#include "GSPIEventCube.hpp"
#include "GSPIEventBin.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_IRF           "GSPIResponse::irf(GEvent&, GSource&, GObservation&)"
#define G_NROI            "GSPIResponse::nroi(GModelSky&, GEnergy&, GTime&, "\
                                                             "GObservation&)"
#define G_EBOUNDS                           "GSPIResponse::ebounds(GEnergy&)"
#define G_SET                 "GSPIResponse::set(const GSPIObservation& obs)"
#define G_LOAD_IRF                       "GSPIResponse::load_irf(GFilename&)"
#define G_LOAD_IRFS                           "GSPIResponse::load_irfs(int&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
//#define G_DEBUG_COMPUTE_IRF                 //!< Debug compute_irf() method

/* __ Constants __________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                       Constructors/destructors                          =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 *
 * Creates an empty INTEGRAL/SPI response.
 ***************************************************************************/
GSPIResponse::GSPIResponse(void) : GResponse()
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Response constructor
 *
 * Constructs an INTEGRAL/SPI response for an observation using a response
 * group filename.
 *
 * This constructor simply stores the file name of a response group, the
 * actual loading will be done using the set() method.
 ***************************************************************************/
GSPIResponse::GSPIResponse(const GFilename& rspname) : GResponse()
{
    // Initialise members
    init_members();

    // Store response group filename
    m_rspname = rspname;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] rsp INTEGRAL/SPI response.
 **************************************************************************/
GSPIResponse::GSPIResponse(const GSPIResponse& rsp) : GResponse(rsp)
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
GSPIResponse::~GSPIResponse(void)
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
 * @param[in] rsp INTEGRAL/SPI response.
 * @return INTEGRAL/SPI response.
 *
 * Assign INTEGRAL/SPI response to this object. The assignment performs
 * a deep copy of all information, hence the original object from which the
 * assignment was made can be destroyed after this operation without any loss
 * of information.
 ***************************************************************************/
GSPIResponse& GSPIResponse::operator=(const GSPIResponse& rsp)
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
 * Clears INTEGRAL/SPI response by resetting all class members to an initial
 * state. Any information that was present before will be lost.
 ***************************************************************************/
void GSPIResponse::clear(void)
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
 * @return Pointer to deep copy of INTEGRAL/SPI response.
 ***************************************************************************/
GSPIResponse* GSPIResponse::clone(void) const
{
    return new GSPIResponse(*this);
}


/***********************************************************************//**
 * @brief Return value of INTEGRAL/SPI instrument response for a photon
 *
 * @param[in] event Observed event.
 * @param[in] photon Incident photon.
 * @param[in] obs Observation.
 * @return Instrument response \f$(cm^2 sr^{-1})\f$
 *
 * @exception GException::invalid_argument
 *            Observation is not a INTEGRAL/SPI observation.
 *
 * Returns the instrument response function for a given observed photon
 * direction as function of the assumed true photon direction. The result
 * is given by
 *
 * \f[
 *    R(p'|p) = 
 * \f]
 *
 * @todo Write down formula
 * @todo Describe in detail how the response is computed.
 ***************************************************************************/
double GSPIResponse::irf(const GEvent&       event,
                         const GPhoton&      photon,
                         const GObservation& obs) const
{
    // Initialise IRF
    double irf = 0.0;

    // Extract INTEGRAL/SPI observation
    const GSPIObservation* spi_obs = dynamic_cast<const GSPIObservation*>(&obs);
    if (spi_obs == NULL) {
        std::string cls = std::string(typeid(&obs).name());
        std::string msg = "Observation of type \""+cls+"\" is not an "
                          "INTEGRAL/SPI observation. Please specify an "
                          "INTEGRAL/SPI observation as argument.";
        throw GException::invalid_argument(G_IRF, msg);
    }

    // Extract INTEGRAL/SPI event cube
    const GSPIEventCube* cube = dynamic_cast<const GSPIEventCube*>(spi_obs->events());
    if (cube == NULL) {
        std::string cls = std::string(typeid(&obs).name());
        std::string msg = "INTEGRAL/SPI observation does not contain a valid "
                          "event cube. Please specify an observation with an "
                          "event cube as argument.";
        throw GException::invalid_argument(G_IRF, msg);
    }

    // Extract INTEGRAL/SPI event bin
    const GSPIEventBin* bin = dynamic_cast<const GSPIEventBin*>(&event);
    if (bin == NULL) {
        std::string cls = std::string(typeid(&event).name());
        std::string msg = "Event of type \""+cls+"\" is  not an INTEGRAL/SPI "
                          "event. Please specify an INTEGRAL/SPI event as "
                          "argument.";
        throw GException::invalid_argument(G_IRF, msg);
    }

    // Continue only if livetime for event is positive
    if (bin->livetime() > 0.0) {

        // Get energy bins of event and photon. Since only the photo peak
        // response is supported so far, both indices need to be identical.
        // Continue only if this is the case.
        if (cube->ebounds().index(photon.energy()) == bin->iebin()) {

            // Get IRF value for photo peak
            irf = irf_value(photon.dir(), *bin, 0);

            // Compile option: Check for NaN/Inf
            #if defined(G_NAN_CHECK)
            if (gammalib::is_notanumber(irf) || gammalib::is_infinite(irf)) {
                std::cout << "*** ERROR: GSPIResponse::irf:";
                std::cout << " NaN/Inf encountered";
                std::cout << " (irf=" << irf;
                std::cout << ")";
                std::cout << std::endl;
            }
            #endif

        } // endif: photon energy in same energy bin as event energy

    } // endif: livetime of event was positive

    // Return IRF value
    return irf;
}


/***********************************************************************//**
 * @brief Return value of INTEGRAL/SPI instrument response for sky direction
 *        and event bin
 *
 * @param[in] srcDir Sky direction.
 * @param[in] bin INTEGRAL/SPI event bin.
 * @param[in] ireg IRF region (0: photo peak).
 * @return Instrument response \f$(cm^2 sr^{-1})\f$
 *
 * Returns the instrument response function for a given sky direction and
 * event bin. The value of the IRF is bilinearly interpolated from the
 * pre-computed IRFs cube that is stored in the class.
 ***************************************************************************/
double GSPIResponse::irf_value(const GSkyDir&      srcDir,
                               const GSPIEventBin& bin,
                               const int&          ireg) const
{
    // Initialise IRF value
    double irf = 0.0;

    // Convert sky direction to zenith angle. Continue only is zenith angle
    // is below maximum zenith angle
    double zenith  = this->zenith(bin.ipt(),  srcDir);
    if (zenith < m_max_zenith) {

        // Convert sky direction to azimuth angle
        double azimuth = this->azimuth(bin.ipt(), srcDir);

        // Compute pixel
        double xpix = (zenith * std::cos(azimuth) - m_wcs_xmin) / m_wcs_xbin;
        double ypix = (zenith * std::sin(azimuth) - m_wcs_ymin) / m_wcs_ybin;

        // Continue only if pixel is within IRF
        if (xpix > 0.0 && xpix < m_wcs_xpix_max &&
            ypix > 0.0 && ypix < m_wcs_ypix_max) {

            // Get number of pixels in X direction
            int nx   = m_irfs.nx();
            int ndet = m_irfs.shape()[0];
            int nreg = m_irfs.shape()[1];

            // Get IRF detector for event
            int idet = irf_detid(bin.dir().detid());
            int ieng = bin.iebin();

            // Get map
            int map = idet + (ireg + ieng * nreg) * ndet;

            // Get 4 nearest neighbours
            int ix_left   = int(xpix);
            int ix_right  = ix_left + 1;
            int iy_top    = int(ypix);
            int iy_bottom = iy_top  + 1;

            // Get weighting factors
            double wgt_right  = xpix - double(ix_left);
            double wgt_left   = 1.0  - wgt_right;
            double wgt_bottom = ypix - double(iy_top);
            double wgt_top    = 1.0  - wgt_bottom;

            // Get indices of 4 nearest neighbours
            int inx_1 = ix_left  + iy_top    * nx;
            int inx_2 = ix_right + iy_top    * nx;
            int inx_3 = ix_left  + iy_bottom * nx;
            int inx_4 = ix_right + iy_bottom * nx;

            // Get weights of 4 nearest neighbours
            double wgt_1 = wgt_left  * wgt_top;
            double wgt_2 = wgt_right * wgt_top;
            double wgt_3 = wgt_left  * wgt_bottom;
            double wgt_4 = wgt_right * wgt_bottom;

            // Compute IRF
            irf  = m_irfs(inx_1, map) * wgt_1;
            irf += m_irfs(inx_2, map) * wgt_2;
            irf += m_irfs(inx_3, map) * wgt_3;
            irf += m_irfs(inx_4, map) * wgt_4;

            // Make sure that IRF does not get negative
            if (irf < 0.0) {
                irf = 0.0;
            }

        } // endif: zenith angle was valid

    } // endif: pixel was within IRF

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
double GSPIResponse::nroi(const GModelSky&    model,
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
 *
 * @todo Implement this method if you need energy dispersion.
 ***************************************************************************/
GEbounds GSPIResponse::ebounds(const GEnergy& obsEnergy) const
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
 * @brief Set response for a specific observation
 *
 * @param[in] obs INTEGRAL/SPI observation.
 * @param[in] energy Line energy
 *
 * Set response for a specific INTEGRAL/SPI observation.
 *
 * If the @p energy argument is set to a positive value, the IRF will be
 * computed for the specified line energy.
 ***************************************************************************/
void GSPIResponse::set(const GSPIObservation& obs, const GEnergy& energy)
{
    // Reset response information
    m_detids.clear();
    m_energies.clear();
    m_ebounds.clear();
    m_irfs.clear();
    m_has_wcs = false;

    // Continue only if the observation contains an event cube
    const GSPIEventCube* cube = dynamic_cast<const GSPIEventCube*>(obs.events());
    if (cube != NULL) {

        // Set requested detector identifiers
        set_detids(cube);

        // Set coordination transformation cache
        set_cache(cube);

        // Load IRFs for photo peak region (region 0)
        load_irfs(0);

        // Get dimension of pre-computed IRF
        int npix = m_irfs.npix();
        int ndet = m_irfs.shape()[0];
        int nreg = m_irfs.shape()[1];
        int neng = cube->naxis(2);

        // Initialise pre-computed IRF
        GSkyMap irfs = m_irfs;
        irfs.nmaps(ndet * nreg * neng);
        irfs.shape(ndet, nreg, neng);
        irfs = 0.0;

        // Pre-compute IRFs for all energy bins
        for (int ieng = 0; ieng < neng; ++ieng) {

            // Get energy boundaries. If line argument is true then use the
            double emin = cube->ebounds().emin(ieng).keV();
            double emax = cube->ebounds().emax(ieng).keV();

            // If line energy is specified then compute a line response for
            // the energy bin that overlaps with the line energy
            double energy_keV = energy.keV();
            if (energy_keV > 0.0) {
                if (energy_keV < emin || energy_keV > emax) {
                    continue;
                }
                emin = energy_keV;
                emax = energy_keV;
            }

            // Pre-compute IRF for this energy bin
            GSkyMap irf = compute_irf(emin, emax);

            // Copy over IRF
            for (int idet = 0; idet < ndet; ++idet) {
                for (int ireg = 0; ireg < nreg; ++ireg) {
                    int map_irf  = idet + ireg * ndet;
                    int map_irfs = idet + (ireg + ieng * nreg) * ndet;
                    for (int i = 0; i < npix; ++i) {
                        irfs(i, map_irfs) = irf(i, map_irf);
                    }
                }
            }

            // Append energy boundary
            m_ebounds.append(GEnergy(emin, "keV"), GEnergy(emax, "keV"));

        } // endfor: looped over all energy bins

        // Replace IRFs
        m_irfs = irfs;

        // Clear energies
        m_energies.clear();

    } // endif: observation contained an event cube

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read SPI response from FITS object
 *
 * @param[in] fits FITS object.
 *
 * Reads the SPI response from FITS object.
 ***************************************************************************/
void GSPIResponse::read(const GFits& fits)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Write SPI response into FITS object
 *
 * @param[in,out] fits FITS object.
 *
 * Writes the SPI response into FITS object.
 ***************************************************************************/
void GSPIResponse::write(GFits& fits) const
{
    // Write IRFs as sky map
    m_irfs.write(fits);

    // Write detector identifiers
    write_detids(fits);

    // Write energies
    write_energies(fits);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load SPI response from file
 *
 * @param[in] filename Response file name.
 *
 * Loads SPI response from response file.
 ***************************************************************************/
void GSPIResponse::load(const GFilename& filename)
{
    // Open FITS file
    GFits fits(filename);

    // Read response
    read(fits);

    // Close FITS file
    fits.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save SPI response into file
 *
 * @param[in] filename Response file name.
 * @param[in] clobber Overwrite existing FITS file?
 *
 * Saves SPI response into response file.
 ***************************************************************************/
void GSPIResponse::save(const GFilename& filename, const bool& clobber) const
{
    // Creat FITS file
    GFits fits;

    // Write response
    write(fits);

    // Save FITS file
    fits.saveto(filename, clobber);

    // Close FITS file
    fits.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print INTEGRAL/SPI response information
 *
 * @param[in] chatter Chattiness.
 * @return String containing INTEGRAL/SPI response information.
 ***************************************************************************/
std::string GSPIResponse::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Extract information
        int nx   = m_irfs.nx();
        int ny   = m_irfs.ny();
        int ndet = (m_irfs.shape().size() > 0) ? m_irfs.shape()[0] : 0;
        int nreg = (m_irfs.shape().size() > 1) ? m_irfs.shape()[1] : 0;
        int neng = (m_irfs.shape().size() > 2) ? m_irfs.shape()[2] : 0;

        // Append header
        result.append("=== GSPIResponse ===");

        // Append information
        result.append("\n"+gammalib::parformat("Response group filename"));
        result.append(m_rspname.url());
        result.append("\n"+gammalib::parformat("Response map pixels"));
        result.append(gammalib::str(nx)+" * "+gammalib::str(ny));
        result.append("\n"+gammalib::parformat("X axis range"));
        result.append("["+gammalib::str(m_wcs_xmin * gammalib::rad2deg));
        result.append(","+gammalib::str(m_wcs_xmax * gammalib::rad2deg));
        result.append("] deg");
        result.append("\n"+gammalib::parformat("Y axis range"));
        result.append("["+gammalib::str(m_wcs_ymin * gammalib::rad2deg));
        result.append(","+gammalib::str(m_wcs_ymax * gammalib::rad2deg));
        result.append("] deg");
        result.append("\n"+gammalib::parformat("Maximum zenith angle"));
        result.append(gammalib::str(m_max_zenith * gammalib::rad2deg)+" deg");
        result.append("\n"+gammalib::parformat("Number of detectors"));
        result.append(gammalib::str(ndet));
        result.append("\n"+gammalib::parformat("Number of regions"));
        result.append(gammalib::str(nreg));
        result.append("\n"+gammalib::parformat("Number of energies"));
        result.append(gammalib::str(neng));
        result.append("\n"+gammalib::parformat("Continuum IRF gamma"));
        result.append(gammalib::str(m_gamma));
        result.append("\n"+gammalib::parformat("Continnum IRF log(E) step"));
        result.append(gammalib::str(m_dlogE));

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
void GSPIResponse::init_members(void)
{
    // Initialise members
    m_rspname.clear();
    m_detids.clear();
    m_energies.clear();
    m_ebounds.clear();
    m_irfs.clear();
    m_dlogE = 0.03;
    m_gamma = 2.0;

    // Initialise cache
    m_spix.clear();
    m_posang.clear();
    m_has_wcs      = false;
    m_wcs_xmin     = 0.0;
    m_wcs_ymin     = 0.0;
    m_wcs_xmax     = 0.0;
    m_wcs_ymax     = 0.0;
    m_wcs_xbin     = 0.0;
    m_wcs_ybin     = 0.0;
    m_wcs_xpix_max = 0.0;
    m_wcs_ypix_max = 0.0;
    m_max_zenith   = 180.0 * gammalib::deg2rad;

    // Set energy scale for response cache to MeV
    //m_irf_cache.energy_scale(GEnergy(1.0, "MeV"));
    //m_nroi_cache.energy_scale(GEnergy(1.0, "MeV"));

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] rsp INTEGRAL/SPI response function.
 ***************************************************************************/
void GSPIResponse::copy_members(const GSPIResponse& rsp)
{
    // Copy members
    m_rspname  = rsp.m_rspname;
    m_detids   = rsp.m_detids;
    m_energies = rsp.m_energies;
    m_ebounds  = rsp.m_ebounds;
    m_irfs     = rsp.m_irfs;
    m_dlogE    = rsp.m_dlogE;
    m_gamma    = rsp.m_gamma;

    // Copy cache
    m_spix         = rsp.m_spix;
    m_posang       = rsp.m_posang;
    m_has_wcs      = rsp.m_has_wcs;
    m_wcs_xmin     = rsp.m_wcs_xmin;
    m_wcs_ymin     = rsp.m_wcs_ymin;
    m_wcs_xmax     = rsp.m_wcs_xmax;
    m_wcs_ymax     = rsp.m_wcs_ymax;
    m_wcs_xbin     = rsp.m_wcs_xbin;
    m_wcs_ybin     = rsp.m_wcs_ybin;
    m_wcs_xpix_max = rsp.m_wcs_xpix_max;
    m_wcs_ypix_max = rsp.m_wcs_ypix_max;
    m_max_zenith   = rsp.m_max_zenith;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GSPIResponse::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Write detector identifiers into FITS object
 *
 * @param[in,out] fits FITS object.
 *
 * Writes the detector identifiers stored into the m_detids member into a
 * FITS object.
 ***************************************************************************/
void GSPIResponse::write_detids(GFits& fits) const
{
    // Get number of detector identifiers
    int num = m_detids.size();

    // Create columns
    GFitsTableShortCol col_detid("DET_ID", num);

    // Fill energy column in units of keV
    for (int i = 0; i < num; ++i) {
        col_detid(i) = m_detids[i];
    }

    // Create energies table
    GFitsBinTable table(num);
    table.append(col_detid);
    table.extname("DETIDS");

    // Append detector identifiers to FITS file
    fits.append(table);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write energies into FITS object
 *
 * @param[in,out] fits FITS object.
 *
 * Writes the energy nodes stored into the m_energies member into a FITS
 * object. The energies are written in units of keV.
 ***************************************************************************/
void GSPIResponse::write_energies(GFits& fits) const
{
    // If there are energy boundaries then the IRF was precomputed and
    // hence we store the energy boundaries in the FITS file
    if (!m_ebounds.is_empty()) {
        m_ebounds.write(fits);
    }

    // ... otherwise we write the energies
    else {

        // Get number of energies
        int num = m_energies.size();

        // Create columns
        GFitsTableDoubleCol col_energy("ENERGY", num);

        // Fill energy column in units of keV
        for (int i = 0; i < num; ++i) {
            col_energy(i) = m_energies[i];
        }
        col_energy.unit("keV");

        // Create energies table
        GFitsBinTable table(num);
        table.append(col_energy);
        table.extname(gammalib::extname_energies);

        // Append energies table to FITS file
        fits.append(table);

    } // endelse: energies written

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load IRF as sky map
 *
 * @param[in] irfname IRF file name.
 *
 * Loads an IRF FITS file as a sky map. The sky map is returned in ARC
 * projection that is a zenith equidistant projection, which is the
 * projection that is natively used to store the IRFs.
 ***************************************************************************/
GSkyMap GSPIResponse::load_irf(const GFilename& irfname) const
{
    // Open IRF FITS file
    GFits fits(irfname);

    // Access IRF FITS image
    GFitsImage* image = fits.image("SPI.-IRF.-RSP");

    // Get image attributes
    int    naxis1 = image->integer("NAXIS1");
    int    naxis2 = image->integer("NAXIS2");
    int    naxis3 = image->integer("NAXIS3");
    int    naxis4 = image->integer("NAXIS4");
    double crval2 = image->real("CRVAL2");
    double crval3 = image->real("CRVAL3");
    double crpix2 = image->real("CRPIX2");
    double crpix3 = image->real("CRPIX3");
    double cdelt2 = image->real("CDELT2");
    double cdelt3 = image->real("CDELT3");

    // Derive image limits. Limits are taken at the pixel centres since
    // we want to use them for bilinear interpolation. This means that
    // we will throw away half a pixel at the edge of the IRFs.
    double wcs_xmin     = (crval2 - (crpix2-1.0) * cdelt2) * gammalib::deg2rad;
    double wcs_ymin     = (crval3 - (crpix3-1.0) * cdelt3) * gammalib::deg2rad;
    double wcs_xbin     = cdelt2 * gammalib::deg2rad;
    double wcs_ybin     = cdelt3 * gammalib::deg2rad;
    double wcs_xmax     = wcs_xmin + double(naxis2-1) * wcs_xbin;
    double wcs_ymax     = wcs_ymin + double(naxis3-1) * wcs_ybin;
    double wcs_xpix_max = double(naxis2-1);
    double wcs_ypix_max = double(naxis3-1);

    // If no image limits exists so far then store them for fast IRF access
    // that does not depend on the actual IRF projection (we just use here
    // the sky map as a convient container)
    if (!m_has_wcs) {
        m_has_wcs      = true;
        m_wcs_xmin     = wcs_xmin;
        m_wcs_ymin     = wcs_ymin;
        m_wcs_xbin     = wcs_xbin;
        m_wcs_ybin     = wcs_ybin;
        m_wcs_xmax     = wcs_xmax;
        m_wcs_ymax     = wcs_ymax;
        m_wcs_xpix_max = wcs_xpix_max;
        m_wcs_ypix_max = wcs_ypix_max;
    }

    // ... otherwise check if the limits are consistent
    else {
        if ((std::abs(wcs_xmin     - m_wcs_xmin)     > 1.0e-6) ||
            (std::abs(wcs_ymin     - m_wcs_ymin)     > 1.0e-6) ||
            (std::abs(wcs_xbin     - m_wcs_xbin)     > 1.0e-6) ||
            (std::abs(wcs_ybin     - m_wcs_ybin)     > 1.0e-6) ||
            (std::abs(wcs_xmax     - m_wcs_xmax)     > 1.0e-6) ||
            (std::abs(wcs_ymax     - m_wcs_ymax)     > 1.0e-6) ||
            (std::abs(wcs_xpix_max - m_wcs_xpix_max) > 1.0e-6) ||
            (std::abs(wcs_ypix_max - m_wcs_ypix_max) > 1.0e-6)) {
            std::string msg = "Inconsistent IRFs encountered in file \""+
                              irfname.url()+"\". Please specify a response "
                              "group where all IRFs have the same definition.";
            throw GException::invalid_value(G_LOAD_IRF, msg);
        }
    }

    // Set maximum zenith angle
    m_max_zenith = (std::abs(m_wcs_xmax) > std::abs(m_wcs_ymax)) ?
                    std::abs(m_wcs_xmax) : std::abs(m_wcs_ymax);

    // Compute sky map attributes
    int nmap = naxis1 * naxis4;

    // Initialise IRF in celestial coordinates. We use an ARC projection as
    // this is the native projection of the IRF files.
    GSkyMap irf("ARC", "CEL", crval2, crval3, cdelt2, cdelt3, naxis2, naxis3, nmap);

    // Fill sky map
    for (int ix = 0; ix < naxis2; ++ix) {
        for (int iy = 0; iy < naxis3; ++iy) {
            for (int idet = 0; idet < naxis1; ++idet) {
                for (int ireg = 0; ireg < naxis4; ++ireg) {
                    int index      = ix + iy * naxis2;
                    int map        = idet + ireg * naxis1;
                    irf(index,map) = image->pixel(idet, ix, iy, ireg);
                }
            }
        }
    }

    // Set sky map shape
    irf.shape(naxis1, naxis4);

    // Close FITS file
    fits.close();

    // Return IRF
    return irf;
}


/***********************************************************************//**
 * @brief Load Instrument Response Functions
 *
 * @param[in]  region IRF region (-1 = load all regions)
 *
 * The method requires that the required detector identifiers were previously
 * setup using the set_detids() method.
 ***************************************************************************/
void GSPIResponse::load_irfs(const int& region)
{
    // Initialise results
    m_energies.clear();
    m_irfs.clear();

    // Determine number of requested detectors. Continue only if there are
    // required detectors.
    int num_det = m_detids.size();
    if (num_det > 0) {

        // Open response group
        GFits rsp(m_rspname);

        // Get response grouping table
        const GFitsTable* grp = gammalib::spi_hdu(rsp, "SPI.-IRF.-RSP-IDX");
        if (grp == NULL) {
            std::string msg = "No response grouping table found in file \""+
                              m_rspname.url()+"\". Please specify a valid "
                              "response grouping file.";
            throw GException::invalid_value(G_LOAD_IRFS, msg);
        }

        // Determine number of response members
        int num_irfs = gammalib::spi_num_hdus(rsp, "SPI.-IRF.-RSP");

        // Setup node array for response energies
        for (int i_irf = 0; i_irf < num_irfs; ++i_irf) {
            m_energies.append((*grp)["ENERGY"]->real(i_irf));
        }

        // Loop over all IRFs
        for (int i_irf = 0; i_irf < num_irfs; ++i_irf) {

            // Get IRF filename
            std::string irfname = m_rspname.path() +
                                  (*grp)["MEMBER_LOCATION"]->string(i_irf);

            // Load IRF
            GSkyMap irf = load_irf(irfname);

            // Determine number of requested regions
            int num_regions = (region == -1) ? irf.shape()[1] : 1;

            // If this is the first IRF then allocate IRFs
            if (i_irf == 0) {

                // Set IRFs from IRF
                m_irfs = irf;

                // Allocate sufficient maps and reshape maps
                m_irfs.nmaps(num_det * num_regions * num_irfs);
                m_irfs.shape(num_det, num_regions, num_irfs);

                // Initialise maps
                m_irfs = 0.0;

            } // endif: this was the first IRF

            // Extract relevant part or IRF
            int nx   = irf.nx();
            int ny   = irf.ny();
            int ndet = irf.shape()[0];
            int nreg = irf.shape()[1];

            // Loop over requested regions
            for (int i_region = 0; i_region < num_regions; ++i_region) {

                // Loop over requested detectors
                for (int i_det = 0; i_det < num_det; ++i_det) {

                    // Set map indices in source IRF and target IRFs
                    int map_irf  = m_detids[i_det] + i_region * ndet;
                    int map_irfs = i_det + (i_region + i_irf * num_regions) * num_det;

                    // Copy IRF pixels
                    for (int i = 0, iy = 0; iy < ny; ++iy) {
                        for (int ix = 0; ix < nx; ++ix, ++i) {
                            m_irfs(i, map_irfs) = irf(i, map_irf);
                        }
                    }

                } // endfor: looped over requested detectors

            } // endfor: looped over requested regions

        } // endfor: looped over all IRFs

        // Close FITS file
        rsp.close();

    } // endif: there were requested detectors

    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute as sky map
 *
 * @param[in] emin Minimum energy (keV).
 * @param[in] emax Maximum energy (keV).
 *
 * Compute IRF for an energy band specified by an energy interval. If the
 * width of the energy interval is zero a line IRF will be computed.
 ***************************************************************************/
GSkyMap GSPIResponse::compute_irf(const double& emin, const double& emax) const
{
    // Get IRF dimension
    int npix = m_irfs.npix();
    int ndet = m_irfs.shape()[0];
    int nreg = m_irfs.shape()[1];

    // Debug: print dimensions
    #if defined(G_DEBUG_COMPUTE_IRF)
    std::cout << "GSPIResponse::compute_irf:" << std::endl;
    std::cout << "- emin = " << emin << " keV" << std::endl;
    std::cout << "- emax = " << emax << " keV" << std::endl;
    std::cout << "- npix = " << npix << std::endl;
    std::cout << "- ndet = " << ndet << std::endl;
    std::cout << "- nreg = " << nreg << std::endl;
    #endif

    // Initialise IRF
    GSkyMap irf = m_irfs;
    irf.nmaps(ndet * nreg);
    irf.shape(ndet, nreg);
    irf = 0.0;

    // Case A: Integrate response over energy band
    if (emax > emin) {

        // Get logarithmic energy limits
        double log_emin   = std::log10(emin);
        double log_emax   = std::log10(emax);
        double log_ewidth = log_emax - log_emin;

        // Set energy weighting factors
        double beta = 1.0 - m_gamma;
        double wgt0 = 1.0 / irf_weight(beta, log_emin, log_emax);

        // Determine number of integration intervals and logarithmic step
        // size
        int num_int = (m_dlogE > DBL_MIN) ? int(log_ewidth / m_dlogE + 0.5) : 1;
        if (num_int < 1) {
            num_int = 1;
        }
        double log_estep = log_ewidth / num_int;

        // Debug: print weighting and internal binning
        #if defined(G_DEBUG_COMPUTE_IRF)
        std::cout << "- log_emin = " << log_emin << std::endl;
        std::cout << "- log_emax = " << log_emax << std::endl;
        std::cout << "- log_ewidth = " << log_ewidth << std::endl;
        std::cout << "- beta = " << beta << std::endl;
        std::cout << "- wgt0 = " << wgt0 << std::endl;
        std::cout << "- num_int = " << num_int << std::endl;
        std::cout << "- log_estep = " << log_estep << std::endl;
        double wgt_check = 0.0;
        #endif

        // Loop over integration intervals
        for (int i_int = 0; i_int < num_int; ++i_int) {

            // Determine integration interval width and midpoint
            double log_emin_bin = log_emin + i_int * log_estep;
            double log_emax_bin = log_emin_bin + log_estep;
            double energy       = std::exp(0.5 * gammalib::ln10 *
                                           (log_emin_bin + log_emax_bin));

            // Get response weighting factor
            double wgt = irf_weight(beta, log_emin_bin, log_emax_bin) * wgt0;

            // Get IRFs bracketing the mean energy and the IRF weights
            m_energies.set_value(energy);
            int irf_low    = m_energies.inx_left();
            int irf_up     = m_energies.inx_right();
            double wgt_low = m_energies.wgt_left()  * wgt;
            double wgt_up  = m_energies.wgt_right() * wgt;

            // Debug: add weight for final check
            #if defined(G_DEBUG_COMPUTE_IRF)
            wgt_check += wgt_low + wgt_up;
            #endif

            // Add contribution to all IRF pixels
            for (int idet = 0; idet < ndet; ++idet) {
                for (int ireg = 0; ireg < nreg; ++ireg) {
                    int map     = idet + ireg * ndet;
                    int map_low = map + irf_low * ndet * nreg;
                    int map_up  = map + irf_up  * ndet * nreg;
                    for (int i = 0; i < npix; ++i) {
                        irf(i, map) += wgt_low * m_irfs(i, map_low) +
                                       wgt_up  * m_irfs(i, map_up);
                    }
                }
            }

        } // endfor: looped over integration intervals

        // Debug: print weighting check
        #if defined(G_DEBUG_COMPUTE_IRF)
        std::cout << "- wgt_check = " << wgt_check << " (should be 1)" << std::endl;
        #endif

    } // endif: integrated response over energy band

    // Case B: compute response for line energy
    else {

        // Get IRFs bracketing the minimum energy and the IRF weights
        m_energies.set_value(emin);
        int irf_low    = m_energies.inx_left();
        int irf_up     = m_energies.inx_right();
        double wgt_low = m_energies.wgt_left();
        double wgt_up  = m_energies.wgt_right();

        // Compute all IRF pixels
        for (int idet = 0; idet < ndet; ++idet) {
            for (int ireg = 0; ireg < nreg; ++ireg) {
                int map     = idet + ireg * ndet;
                int map_low = map + irf_low * ndet * nreg;
                int map_up  = map + irf_up  * ndet * nreg;
                for (int i = 0; i < npix; ++i) {
                    irf(i, map) = wgt_low * m_irfs(i, map_low) +
                                  wgt_up  * m_irfs(i, map_up);
                }
            }
        }

        // Debug: print weighting check
        #if defined(G_DEBUG_COMPUTE_IRF)
        std::cout << "- wgt_low = " << wgt_low << std::endl;
        std::cout << "- wgt_up = " << wgt_up << std::endl;
        std::cout << "- wgt_low+wgt_up = " << wgt_low+wgt_up;
        std::cout << " (should be 1)" << std::endl;
        #endif

    } // endelse: compute response for line energy

    // Return IRF
    return irf;
}


/***********************************************************************//**
 * @brief Set vector of detector identifiers used by the observation
 *
 * @param[in] cube INTEGRAL/SPI event cube.
 *
 * Sets the vector of detector identifiers that is used by the observation.
 * The method scans the detector identifiers for all pointings and builds
 * a vector of unique detector identifiers in the order they were encountered
 * in the event cube.
 ***************************************************************************/
void GSPIResponse::set_detids(const GSPIEventCube* cube)
{
    // Clear vector of detector identifiers
    m_detids.clear();

    // Extract relevant event cube dimensions
    int npt  = cube->naxis(0);
    int ndet = cube->naxis(1);

    // Loop over all pointings and detectors
    for (int ipt = 0; ipt < npt; ++ipt) {
        for (int idet = 0; idet < ndet; ++idet) {

            // Get IRF detector identifier
            int detid = irf_detid(cube->dir(ipt, idet).detid());

            // Push IRF detector identifier on vector if it does not yet exist
            std::vector<int>::iterator it = find(m_detids.begin(), m_detids.end(), detid);
            if (it == m_detids.end()) {
                m_detids.push_back(detid);
            }

        } // endfor: looped over detectors
    } // endfor: looped over pointings

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set computation cache
 *
 * @param[in] cube INTEGRAL/SPI event cube.
 *
 * Setup of two vectors for fast coordinate transformation into the
 * instrument system. The first vector m_spix stores the SPI pointing
 * direction (the X direction) while the second vector stores the position
 * angle in celestial coordinates of the SPI Y direction.
 ***************************************************************************/
void GSPIResponse::set_cache(const GSPIEventCube* cube)
{
    // Extract relevant event cube dimensions
    int npt = cube->naxis(0);

    // Clear vectors of SPI X direction and position angles
    m_spix.clear();
    m_spix.reserve(npt);
    m_posang.clear();
    m_posang.reserve(npt);

    // Loop over all pointings
    for (int ipt = 0; ipt < npt; ++ipt) {

        // Compute position angle
        double posang = cube->spi_x(ipt).posang(cube->spi_z(ipt)) + gammalib::pihalf;

        // Store angles
        m_spix.push_back(cube->spi_x(ipt));
        m_posang.push_back(posang);

    } // endfor: looped over pointings

    // Return
    return;
}


/***********************************************************************//**
 * @brief Convert detector identifier into IRF detector identifier
 *
 * @param[in] detid SPI event detector identifier.
 * @return IRF detector identifier
 *
 * Converts a SPI event detector identifier into an IRF detector identifier.
 * TODO: Describe how this is done and why
 ***************************************************************************/
int GSPIResponse::irf_detid(const int& detid) const
{
    // Initialise irf detector identifier
    int irf_detid = detid;

    // Put detector identifier in the relevant range
    if (irf_detid >= 123) {
        irf_detid -= 123;
    }
    else if (irf_detid >= 104) {
        irf_detid -= 104;
    }
    else if (irf_detid >= 85) {
        irf_detid -= 85;
    }

    // Return irf detector identifier
    return irf_detid;
}


/***********************************************************************//**
 * @brief Compute weight of logarithmic energy bin
 *
 * @param[in] beta 1-gamma.
 * @param[in] emin Logarithmic minimum energy.
 * @param[in] emax Logarithmic maximum energy.
 ***************************************************************************/
double GSPIResponse::irf_weight(const double& beta,
                                const double& emin,
                                const double& emax) const
{
    // Initialise weight
    double weight;

    // Assign weight
    if (std::abs(beta) < DBL_MIN) {
        weight = (emax - emin);
    }
    else {
        weight = std::exp(beta * gammalib::ln10 * emax) -
                 std::exp(beta * gammalib::ln10 * emin);
    }

    // Return weight
    return weight;
}
