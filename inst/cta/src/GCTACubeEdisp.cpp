/***************************************************************************
 *      GCTACubeEdisp.cpp - CTA cube analysis energy dispersion class      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016 by Michael Mayer                                    *
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
 * @file GCTACubeEdisp.cpp
 * @brief CTA cube analysis energy dispersion class implementation
 * @author Michael Mayer
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GMath.hpp"
#include "GLog.hpp"
#include "GFits.hpp"
#include "GFitsImage.hpp"
#include "GFitsTable.hpp"
#include "GObservations.hpp"
#include "GCTACubeEdisp.hpp"
#include "GCTAObservation.hpp"
#include "GCTAResponseIrf.hpp"
#include "GCTAEventList.hpp"
#include "GCTAEventCube.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_CONSTRUCTOR1         "GCTACubeEdisp(GCTAEventCube&, double&, int&)"
#define G_CONSTRUCTOR2  "GCTACubeEdisp(std::string&, std::string&, double&, "\
                         "double&, double&, double&, int&, int&, GEbounds&, "\
                                                             "double&, int&)"
#define G_EBOUNDS                          "GCTACubeEdisp::ebounds(GEnergy&)"
#define G_FILL_CUBE   "GCTACubeEdisp::fill_cube(GCTAObservation&, GSkyMap*, "\
                                                                     "GLog*)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Constants __________________________________________________________ */
const GEnergy g_energy_margin(1.0e-12, "TeV");


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 *
 * Constructs an empty energy dispersion cube.
 ***************************************************************************/
GCTACubeEdisp::GCTACubeEdisp(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] cube Energy dispersion cube.
 *
 * Constructs an energy dispersion cube by copying another energy dispersion
 * cube.
 ***************************************************************************/
GCTACubeEdisp::GCTACubeEdisp(const GCTACubeEdisp& cube)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(cube);

    // Return
    return;
}


/***********************************************************************//**
 * @brief File constructor
 *
 * @param[in] filename Energy dispersion cube filename.
 *
 * Constructs an energy dispersion cube by loading the energy dispersion
 * information from an energy dispersion cube FITS file. See the load()
 * method for details about the format of the FITS file.
 ***************************************************************************/
GCTACubeEdisp::GCTACubeEdisp(const GFilename& filename)
{
    // Initialise class members
    init_members();

    // Load Edisp cube from file
    load(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Event cube constructor
 *
 * @param[in] cube Event cube.
 * @param[in] mmax Maximum energy migration (degrees, >0.0).
 * @param[in] nmbins Number of migration bins (2 ...).
 *
 * @exception GException::invalid_argument
 *            Maximum energy migration or number of migration bins invalid.
 *
 * Construct an energy dispersion cube with all elements set to zero using
 * the same binning and sky projection that is used for an event cube.
 ***************************************************************************/
GCTACubeEdisp::GCTACubeEdisp(const GCTAEventCube& cube, const double& mmax,
                             const int& nmbins)
{
    // Throw an exception if the number of migra bins is invalid or the
    // maximum migration is not positive
    if (nmbins < 2) {
        std::string msg = "Number "+gammalib::str(nmbins)+" of migration "
                          "bins is smaller than 2. Please request at least "
                          "2 migration bins.";
        throw GException::invalid_argument(G_CONSTRUCTOR1, msg);
    }
    if (mmax <= 0.0) {
        std::string msg = "Maximum migration "+gammalib::str(mmax)+" is not "
                          "positive. Please specify a positive maximum "
                          "energy migration.";
        throw GException::invalid_argument(G_CONSTRUCTOR1, msg);
    }

    // Initialise class members
    init_members();

    // Use the energy bin boundaries of the event cube as the node energies
    // of the energy dispersion cube.
    m_energies.clear();
    int nbins = cube.ebounds().size();
    if (nbins > 0) {
        m_energies.append(cube.ebounds().emin(0));
        for (int i = 0; i < nbins; ++i) {
            m_energies.append(cube.ebounds().emax(i));
        }
    }

    // Set energy node array used for interpolation
    set_eng_axis();

    // Set the migration nodes
    m_migras.clear();
    double binsize = mmax / double(nmbins);
    for (int i = 0; i < nmbins; ++i) {
        double migra = binsize * (double(i) + 0.5); // avoid central singularity
        m_migras.append(migra);
    }

    // Compute number of sky maps
    int nmaps = m_energies.size() * m_migras.size();

    // Set energy dispersion cube to event cube
    m_cube = cube.counts();

    // Set appropriate number of skymaps
    m_cube.nmaps(nmaps);

    // Set cube shape
    m_cube.shape(m_migras.size(), m_energies.size());

    // Set all energy dispersion cube pixels to zero as we want to have
    // a clean map upon construction
    m_cube = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Energy dispersion cube constructor
 *
 * @param[in] wcs      World Coordinate System.
 * @param[in] coords   Coordinate System (CEL or GAL).
 * @param[in] x        X coordinate of sky map centre (deg).
 * @param[in] y        Y coordinate of sky map centre (deg).
 * @param[in] dx       Pixel size in x direction at centre (deg/pixel).
 * @param[in] dy       Pixel size in y direction at centre (deg/pixel).
 * @param[in] nx       Number of pixels in x direction.
 * @param[in] ny       Number of pixels in y direction.
 * @param[in] energies True energies.
 * @param[in] mmax     Maximum energy migration (degrees, >0.0).
 * @param[in] nmbins   Number of migration bins (2 ...).
 *
 * @exception GException::invalid_argument
 *            Maximum energy migration or number of migration bins invalid.
 *
 * Constructs an energy dispersion cube by specifying the sky map grid and
 * the energies.
 ***************************************************************************/
GCTACubeEdisp::GCTACubeEdisp(const std::string&   wcs,
                             const std::string&   coords,
                             const double&        x,
                             const double&        y,
                             const double&        dx,
                             const double&        dy,
                             const int&           nx,
                             const int&           ny,
                             const GEnergies&     energies,
                             const double&        mmax,
                             const int&           nmbins)
{
    // Throw an exception if the number of migra bins is invalid or the
    // maximum migration is not positive
    if (nmbins < 2) {
        std::string msg = "Number "+gammalib::str(nmbins)+" of migration "
                          "bins is smaller than 2. Please request at least "
                          "2 migration bins.";
        throw GException::invalid_argument(G_CONSTRUCTOR2, msg);
    }
    if (mmax <= 0.0) {
        std::string msg = "Maximum migration "+gammalib::str(mmax)+" is not "
                          "positive. Please specify a positive maximum "
                          "energy migration.";
        throw GException::invalid_argument(G_CONSTRUCTOR2, msg);
    }

    // Initialise class members
    init_members();

    // Store true energies
    m_energies = energies;

    // Set energy node array used for interpolation
    set_eng_axis();

    // Set migration node array
    m_migras.clear();
    double binsize = mmax / double(nmbins);
    for (int i = 0; i < nmbins; ++i) {
        double migra = binsize * (double(i) + 0.5); // avoid central singularity
        m_migras.append(migra);
    }

    // Compute number of sky maps
    int nmaps = m_energies.size() * m_migras.size();

    // Create sky map
    m_cube = GSkyMap(wcs, coords, x, y, dx, dy, nx, ny, nmaps);

    // Set cube shape
    m_cube.shape(m_migras.size(), m_energies.size());

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 *
 * Destructs energy dispersion cube.
 ***************************************************************************/
GCTACubeEdisp::~GCTACubeEdisp(void)
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
 * @param[in] cube Energy dispersion cube.
 * @return Energy dispersion cube.
 *
 * Assigns energy dispersion cube.
 ***************************************************************************/
GCTACubeEdisp& GCTACubeEdisp::operator=(const GCTACubeEdisp& cube)
{
    // Execute only if object is not identical
    if (this != &cube) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(cube);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/***********************************************************************//**
 * @brief Return energy dispersion (in units of MeV^-1)
 *
 * @param[in] dir Coordinate of the true photon position.
 * @param[in] migra Energy migration (reconstructed over true energy).
 * @param[in] srcEng True photon energy.
 * @return Energy dispersion (in units of MeV^-1)
 *
 * Returns the energy dispersion for a given energy migration (reconstructed
 * over true energy), true photon energy, and sky direction in units of
 * MeV^-1.
 ***************************************************************************/
double GCTACubeEdisp::operator()(const GSkyDir& dir,
                                 const double&  migra,
                                 const GEnergy& srcEng) const
{
    // Update indices and weighting factors for interpolation
    update(migra, srcEng.log10TeV());

    // Perform bi-linear interpolation
    double edisp = m_wgt1 * m_cube(dir, m_inx1) +
                   m_wgt2 * m_cube(dir, m_inx2) +
                   m_wgt3 * m_cube(dir, m_inx3) +
                   m_wgt4 * m_cube(dir, m_inx4);

    // Make sure that energy dispersion does not become negative
    if (edisp < 0.0) {
        edisp = 0.0;
    }

    // Return energy dispersion
    return edisp;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear energy dispersion cube
 *
 * Clears energy dispersion cube.
 ***************************************************************************/
void GCTACubeEdisp::clear(void)
{
    // Free class members
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone energy dispersion cube
 *
 * @return Pointer to deep copy of energy dispersion cube.
 ***************************************************************************/
GCTACubeEdisp* GCTACubeEdisp::clone(void) const
{
    return new GCTACubeEdisp(*this);
}


/***********************************************************************//**
 * @brief Set energy dispersion cube for one CTA observation
 *
 * @param[in] obs CTA observation.
 *
 * Sets the energy dispersion cube for one CTA observation.
 ***************************************************************************/
void GCTACubeEdisp::set(const GCTAObservation& obs)
{
    // Clear energy dispersion cube
    clear_cube();

    // Fill energy dispersion cube
    fill_cube(obs);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Fill energy dispersion cube from observation container
 *
 * @param[in] obs Observation container.
 * @param[in] log Pointer towards logger.
 *
 * Sets the energy dispersion cube from all CTA observations in the
 * observation container.
 ***************************************************************************/
void GCTACubeEdisp::fill(const GObservations& obs, GLog* log)
{
    // Clear energy dispersion cube
    clear_cube();

    // Initialise skymap for exposure weight accumulation
    GSkyMap exposure(m_cube);

    // Loop over all observations in container
    for (int i = 0; i < obs.size(); ++i) {

        // Get observation and continue only if it is a CTA observation
        const GCTAObservation* cta = dynamic_cast<const GCTAObservation*>
                                     (obs[i]);

        // Skip observation if it's not CTA
        if (cta == NULL) {
            if (log != NULL) {
                *log << "Skipping ";
                *log << obs[i]->instrument();
                *log << " observation ";
                *log << "\"" << obs[i]->name() << "\"";
                *log << " (id=" << obs[i]->id() << ")" << std::endl;
            }
            continue;
        }

        // Skip observation if we have a binned observation
        if (cta->eventtype() == "CountsCube") {
            if (log != NULL) {
                *log << "Skipping binned ";
                *log << cta->instrument();
                *log << " observation ";
                *log << "\"" << cta->name() << "\"";
                *log << " (id=" << cta->id() << ")" << std::endl;
            }
            continue;
        }

        // Fill exposure cube cube
        fill_cube(*cta, &exposure, log);

    } // endfor: looped over observations

    // Compute mean energy dispersion cube by dividing through the weights
    for (int pixel = 0; pixel < m_cube.npix(); ++pixel) {
        for (int iebin = 0; iebin < m_energies.size(); ++iebin) {
            if (exposure(pixel, iebin) > 0.0) {
                double norm = 1.0 / exposure(pixel, iebin);
                for (int imigra = 0; imigra < m_migras.size(); ++imigra) {
                    int imap = offset(imigra, iebin);
                    m_cube(pixel, imap) *= norm;
                }
            }
            else {
                for (int imigra = 0; imigra < m_migras.size(); ++imigra) {
                    int imap = offset(imigra, iebin);
                    m_cube(pixel, imap) = 0.0;
                }
            }
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read energy dispersion cube from FITS file
 *
 * @param[in] fits FITS file.
 *
 * Reads the energy dispersion cube from a FITS file. The energy dispersion
 * cube values are expected in the primary extension of the FITS file, the
 * true energy boundaries are expected in the "ENERGIES" extension, and the
 * migration values are expected in the "MIGRAS" extension.
 ***************************************************************************/
void GCTACubeEdisp::read(const GFits& fits)
{
    // Clear object
    clear();

    // Get HDUs
    const GFitsImage& hdu_edispcube = *fits.image("Primary");
    const GFitsTable& hdu_energies  = *fits.table("ENERGIES");
    const GFitsTable& hdu_migras    = *fits.table("MIGRAS");

    // Read energy dispersion cube
    m_cube.read(hdu_edispcube);

    // Read true energy nodes
    m_energies.read(hdu_energies);

    // Read migration nodes
    m_migras.read(hdu_migras);

    // Set energy node array used for interpolation
    set_eng_axis();

    // Compute true energy boundaries
    compute_ebounds();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write energy dispersion cube into FITS file.
 *
 * @param[in] fits FITS file.
 *
 * Write the energy dispersion cube into a FITS file. The energy dispersion
 * cube values are written into the primary extension of the FITS file, the
 * true energy boundaries are written into the "EBOUNDS" extension, and the
 * migration values are written into the "MIGRAS" extension.
 ***************************************************************************/
void GCTACubeEdisp::write(GFits& fits) const
{
    // Remove extensions from FITS file
    if (fits.contains("MIGRAS")) {
        fits.remove("MIGRAS");
    }
    if (fits.contains("ENERGIES")) {
        fits.remove("ENERGIES");
    }
    if (fits.contains("Primary")) {
        fits.remove("Primary");
    }

    // Write cube
    m_cube.write(fits);

    // Write energy boundaries
    m_energies.write(fits, "ENERGIES");

    // Write migration nodes
    m_migras.write(fits, "MIGRAS");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load energy dispersion cube from FITS file
 *
 * @param[in] filename FITS file name.
 *
 * Loads the energy dispersion cube from a FITS file. See the read() method
 * for information about the expected FITS file structure.
 ***************************************************************************/
void GCTACubeEdisp::load(const GFilename& filename)
{
    // Open FITS file
    GFits fits(filename);

    // Read Edisp cube
    read(fits);

    // Close Edisp file
    fits.close();

    // Store filename
    m_filename = filename;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save energy dispersion cube into FITS file
 *
 * @param[in] filename FITS file name.
 * @param[in] clobber Overwrite existing file?
 *
 * Save the energy dispersion cube into a FITS file.
 ***************************************************************************/
void GCTACubeEdisp::save(const GFilename& filename, const bool& clobber) const
{
    // Create FITS file
    GFits fits;

    // Write PSF cube
    write(fits);

    // Save FITS file
    fits.saveto(filename, clobber);

    // Store filename
    m_filename = filename;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return boundaries in true energy
 *
 * @param[in] obsEng Observed photon energy.
 * @return Boundaries in true energy.
 *
 * @exception GException::invalid_value
 *            No energies defined for energy dispersion cube
 *
 * Returns the boundaries in true photon energies that enclose the energy
 * disperson for a given observed photon energy @p obsEng.
 ***************************************************************************/
GEbounds GCTACubeEdisp::ebounds(const GEnergy& obsEng) const
{
    // Throw an exception if there are no energies
    if (m_energies.is_empty()) {
        std::string msg = "No energies defined for energy dispersion cube. "
                          "Please define the energies before calling the "
                          "method.";
        throw GException::invalid_value(G_EBOUNDS, msg);
    }

	// If ebounds were not computed, before, compute them now
	if (m_ebounds.empty()) {
		compute_ebounds();
	}

    // Return the index of the energy boundaries that are just below the
    // observed energy. As the energy dispersion decreases with increasing
    // energy we assure in that way that we integrate over a sufficiently
    // large true energy range.
    int index = 0;
    if (obsEng > m_energies[0]) {
        for (int i = 1; i < m_energies.size(); ++i) {
            if (obsEng < m_energies[i]) {
                break;
            }
            index++;
        }
        if (index >= m_energies.size()) {
            index = m_energies.size();
        }
    }

	// Return true energy boundaries
	return (m_ebounds[index]);
}


/***********************************************************************//**
 * @brief Print energy dispersion cube information
 *
 * @param[in] chatter Chattiness.
 * @return String containing energy dispersion cube information.
 ***************************************************************************/
std::string GCTACubeEdisp::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCTACubeEdisp ===");

        // Append information
        result.append("\n"+gammalib::parformat("Filename")+m_filename);

        // Append energy intervals
        if (m_ebounds.size() > 0) {
            result.append("\n"+m_energies.print(chatter));
        }
        else {
            result.append("\n"+gammalib::parformat("Energies") +
                          "Not defined");
        }

        // Append number of migra bins
        result.append("\n"+gammalib::parformat("Number of migra bins") +
                      gammalib::str(m_migras.size()));

        // Append migra range
        result.append("\n"+gammalib::parformat("Migra range"));
        if (m_migras.size() > 0) {
            result.append(gammalib::str(m_migras[0]));
            result.append(" - ");
            result.append(gammalib::str(m_migras[m_migras.size()-1]));
        }
        else {
            result.append("not defined");
        }

        // Append skymap definition
        result.append("\n"+m_cube.print(chatter));

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
void GCTACubeEdisp::init_members(void)
{
    // Initialise members
    m_filename.clear();
    m_cube.clear();
    m_energies.clear();
    m_elogmeans.clear();
    m_migras.clear();

    // Initialise cache
    m_inx1 = 0;
    m_inx2 = 0;
    m_inx3 = 0;
    m_inx4 = 0;
    m_wgt1 = 0.0;
    m_wgt2 = 0.0;
    m_wgt3 = 0.0;
    m_wgt4 = 0.0;
    m_ebounds.clear();
   
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] cube Energy dispersion cube.
 ***************************************************************************/
void GCTACubeEdisp::copy_members(const GCTACubeEdisp& cube)
{
    // Copy members
    m_filename  = cube.m_filename;
    m_cube      = cube.m_cube;
    m_energies  = cube.m_energies;
    m_elogmeans = cube.m_elogmeans;
    m_migras    = cube.m_migras;

    // Copy cache
    m_inx1    = cube.m_inx1;
    m_inx2    = cube.m_inx2;
    m_inx3    = cube.m_inx3;
    m_inx4    = cube.m_inx4;
    m_wgt1    = cube.m_wgt1;
    m_wgt2    = cube.m_wgt2;
    m_wgt3    = cube.m_wgt3;
    m_wgt4    = cube.m_wgt4;
    m_ebounds = cube.m_ebounds;

    // Return
    return;
}

/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTACubeEdisp::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Clear all pixels in the energy dispersion cube
 ***************************************************************************/
void GCTACubeEdisp::clear_cube(void)
{
    // Loop over all maps
    for (int imap = 0; imap < m_cube.nmaps(); ++imap) {

        // Loop over all pixels in sky map
        for (int pixel = 0; pixel < m_cube.npix(); ++pixel) {

            // Reset cube value to zero
            m_cube(pixel, imap) = 0.0;

        } // endfor: looped over all pixels

    } // endfor: looped over maps

    // Return
    return;
}


/***********************************************************************//**
 * @brief Fill energy dispersion cube from observation container
 *
 * @param[in] obs Observation.
 * @param[in] exposure Pointer towards exposure map.
 * @param[in] log Pointer towards logger.
 *
 * @exception GException::invalid_value
 *            No RoI or response found in CTA observation.
 ***************************************************************************/
void GCTACubeEdisp::fill_cube(const GCTAObservation& obs,
                              GSkyMap*               exposure,
                              GLog*                  log)
{
    // Only continue if we have an event list
    if (obs.eventtype() == "EventList") {

        // Extract pointing direction, energy boundaries and ROI from
        // observation
        GSkyDir  pnt         = obs.pointing().dir();
        GEbounds obs_ebounds = obs.ebounds();
        GCTARoi  roi         = obs.roi();

        // Check for RoI sanity
        if (!roi.is_valid()) {
            std::string msg = "No RoI information found in "+obs.instrument()+
                              " observation \""+obs.name()+"\". Run ctselect "
                              "to specify an RoI for this observation.";
            throw GException::invalid_value(G_FILL_CUBE, msg);
        }

        // Extract response from observation
        const GCTAResponseIrf* rsp = dynamic_cast<const GCTAResponseIrf*>
                                     (obs.response());
        if (rsp == NULL) {
            std::string msg = "No valid instrument response function found in "+
                              obs.instrument()+" observation \""+obs.name()+
                              "\". Please specify the instrument response "
                              "function for this observation.";
            throw GException::invalid_value(G_FILL_CUBE, msg);
        }

        // Get energy dispersion component
        const GCTAEdisp* edisp = rsp->edisp();
        if (edisp == NULL) {
            std::string msg = "No energy dispersion rcomponent found in the "
                              "the instrument response of "+
                              obs.instrument()+" observation \""+obs.name()+
                              "\". Please provide an instrument response that "
                              "comprises an energy dispersion component.";
            throw GException::invalid_value(G_FILL_CUBE, msg);
        }

        // Announce observation usage
        if (log != NULL) {
            *log << "Including ";
            *log << obs.instrument();
            *log << " observation \"" << obs.name();
            *log << "\" (id=" << obs.id() << ")";
            *log << " in energy dispersion cube computation." << std::endl;
        }

        // Loop over all pixels in sky map
        for (int pixel = 0; pixel < m_cube.npix(); ++pixel) {

            // Get pixel sky direction
            GSkyDir dir = m_cube.inx2dir(pixel);

            // Skip pixel if it is outside the RoI
            if (roi.centre().dir().dist_deg(dir) > roi.radius()) {
                continue;
            }
            
            // Compute theta angle with respect to pointing direction in radians
            double theta = pnt.dist(dir);

            // Loop over all energy bins
            for (int iebin = 0; iebin < m_energies.size(); ++iebin) {

                // Skip energy dispersion cube energy if the energy is
                // outside the observation energy. The energy dispersion cube
                // energies are true energies while the observation energy
                // boundaries are reconstructed energies, hence this is only
                // an approximation, but probably the only we can really do.
                if (!obs_ebounds.contains(m_energies[iebin])) {
                    continue;
                }

                // Get log10 of true energy in TeV
                double logEsrc = m_energies[iebin].log10TeV();

                // Compute exposure weight
                double weight = rsp->aeff(theta, 0.0, 0.0, 0.0, logEsrc) *
                                obs.livetime();

                // If available, accumulate weights
                if (exposure != NULL) {
                    (*exposure)(pixel, iebin) += weight;
                }

                // Loop over delta values
                for (int imigra = 0; imigra < m_migras.size(); ++imigra) {

                    // Compute energy migration fraction
                    double migra = m_migras[imigra];

                    // Skip migra bin if zero
                    if (migra <= 0.0) {
                        continue;
                    }

                    // Compute log10 of reconstructed energy in TeV
                    double  logEobs = std::log10(migra) + logEsrc;
                    GEnergy Eobs(std::pow(10.0, logEobs), "TeV");

                    // Set map index
                    int imap = offset(imigra, iebin);

                    // Add energy dispersion cube value
                    m_cube(pixel, imap) += rsp->edisp(Eobs,
                                                      theta, 0.0,
                                                      0.0, 0.0,
                                                      logEsrc) * weight;

                } // endfor: looped over migration bins
                
            } // endfor: looped over energy bins

        } // endfor: looped over all pixels

    } // endif: observation contained an event list

    // Return
    return;
}


/***********************************************************************//**
 * @brief Update energy dispersion cube indices and weights cache
 *
 * @param[in] migra Reconstructed over true energy.
 * @param[in] logEsrc Log10 of true photon energy in TeV. 
 *
 * Updates the energy dispersion cube indices, stored in the members m_inx1,
 * m_inx2, m_inx3, and m_inx4, and weights cache, stored in m_wgt1, m_wgt2,
 * m_wgt3, and m_wgt4.
 ***************************************************************************/
void GCTACubeEdisp::update(const double& migra, const double& logEsrc) const
{
    // Set node array for energy interpolation
    m_elogmeans.set_value(logEsrc);
   
    // Set node array for delta interpolation
	m_migras.set_value(migra);

    // Set indices for bi-linear interpolation
    m_inx1 = offset(m_migras.inx_left(),  m_elogmeans.inx_left());
    m_inx2 = offset(m_migras.inx_left(),  m_elogmeans.inx_right());
    m_inx3 = offset(m_migras.inx_right(), m_elogmeans.inx_left());
    m_inx4 = offset(m_migras.inx_right(), m_elogmeans.inx_right());

    // Set weighting factors for bi-linear interpolation
    m_wgt1 = m_migras.wgt_left()  * m_elogmeans.wgt_left();
    m_wgt2 = m_migras.wgt_left()  * m_elogmeans.wgt_right();
    m_wgt3 = m_migras.wgt_right() * m_elogmeans.wgt_left();
    m_wgt4 = m_migras.wgt_right() * m_elogmeans.wgt_right();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set nodes for interpolation in true energy
 ***************************************************************************/
void GCTACubeEdisp::set_eng_axis(void)
{
    // Clear nodes
    m_elogmeans.clear();

    // Set nodes
    for (int i = 0; i < m_energies.size(); ++i) {
        m_elogmeans.append(m_energies[i].log10TeV());
    }
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute true energy boundary vector
 *
 * Computes for all energies of the energy dispersion cube the boundaries in
 * true energy that encompass non-negligible migration matrix elements.
 * Values >= 1.0e-12 are considered as non-negligible. In case that no matrix
 * elements are found for a given energy, the interval of true energies will
 * be set to [1 TeV, 1 TeV] (i.e. an empty interval).
 ***************************************************************************/
void GCTACubeEdisp::compute_ebounds() const
{
	// Set epsilon
	const double eps = 1.0e-12;

    // Clear energy boundary vector
    m_ebounds.clear();

    // Loop over all energies
	for (int i = 0; i < m_energies.size(); i++) {

		// Initialise results
		double migra_min = 0.0;
		double migra_max = 0.0;
		bool   minFound  = false;
		bool   maxFound  = false;

		// Compute map of sky ranges to be considered
		int mapmin = m_migras.size() * i;
		int mapmax = m_migras.size() * (i + 1);

		// Loop over sky map entries from lower migra
		for (int imap = mapmin; imap < mapmax; ++imap) {

			// Get maximum energy dispersion term for all sky pixels
            double edisp = 0.0;
            for (int pixel = 0; pixel < m_cube.npix(); ++pixel) {
                double value = m_cube(pixel, imap);
                if (value > edisp) {
                    edisp = value;
                }
            }

			// Find first non-negligible energy dispersion term
			if (edisp >= eps) {
				minFound  = true;
				migra_min = m_migras[imap - mapmin];
				break;
			}

		} // endfor: loop over maps

		// Loop over sky map entries from high migra
		for (int imap = mapmax-1; imap >= mapmin; --imap) {

			// Get maximum energy dispersion term for all sky pixels
            double edisp = 0.0;
            for (int pixel = 0; pixel < m_cube.npix(); ++pixel) {
                double value = m_cube(pixel, imap);
                if (value > edisp) {
                    edisp = value;
                }
            }

			// Find first non-negligible energy dispersion term
			if (edisp >= eps) {
				maxFound  = true;
				migra_max = m_migras[imap - mapmin];
				break;
			}

		} // endfor: loop over maps

		// Initialise energy boundaries
		GEnergy emin;
		GEnergy emax;

		// Compute energy boundaries if they were found and if they are
        // valid
		if (minFound && maxFound && migra_min > 0.0 && migra_max > 0.0 &&
            migra_max > migra_min) {
			emin = m_energies[i] * migra_min;
			emax = m_energies[i] * migra_max;
		}
    
        // ... otherwise we set the interval to a zero interval for safety
		else {
			emin.log10TeV(0.0);
			emax.log10TeV(0.0);
		}

		// Append energy boundaries
		m_ebounds.push_back(GEbounds(emin, emax));

	} // endfor: looped over all energies

    // Return
    return;
}
