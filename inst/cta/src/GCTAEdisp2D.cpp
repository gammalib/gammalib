/***************************************************************************
 *         GCTAEdisp2D.cpp - CTA 2D energy dispersion class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2015 by Juergen Knoedlseder                         *
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
 * @file GCTAEdisp2D.cpp
 * @brief CTA 2D energy dispersion class implementation
 * @author Florent Forest
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include <vector>
#include "GTools.hpp"
#include "GFits.hpp"
#include "GFitsBinTable.hpp"
#include "GFitsTable.hpp"
#include "GCTAEdisp2D.hpp"
#include "GMath.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Constants __________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GCTAEdisp2D::GCTAEdisp2D(void) : GCTAEdisp()
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief File constructor
 *
 * @param[in] filename FITS file name.
 *
 * Construct instance by loading the energy dispersion information from a FITS
 * file.
 ***************************************************************************/
GCTAEdisp2D::GCTAEdisp2D(const std::string& filename) : GCTAEdisp()
{
    // Initialise class members
    init_members();

    // Load energy dispersion from file
    load(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] edisp Energy dispersion.
 ***************************************************************************/
GCTAEdisp2D::GCTAEdisp2D(const GCTAEdisp2D& edisp) : GCTAEdisp(edisp)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(edisp);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCTAEdisp2D::~GCTAEdisp2D(void)
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
 * @param[in] edisp Energy dispersion.
 * @return Energy dispersion.
 ***************************************************************************/
GCTAEdisp2D& GCTAEdisp2D::operator= (const GCTAEdisp2D& edisp)
{
    // Execute only if object is not identical
    if (this != &edisp) {

        // Copy base class members
        this->GCTAEdisp::operator=(edisp);

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(edisp);

    } // endif: object was not identical

    // Return this object
    return *this;
}

/***********************************************************************//**
 * @brief Return energy dispersion in units of s^-1 MeV^-1
 *
 * @param[in] logEobs Log10 of the measured energy (TeV).
 * @param[in] logEsrc Log10 of the true photon energy (TeV).
 * @param[in] theta Offset angle in camera system (rad). Defaults to 0.0.
 * @param[in] phi Azimuth angle in camera system (rad). Not used in this method.
 * @param[in] zenith Zenith angle in Earth system (rad). Not used in this method.
 * @param[in] azimuth Azimuth angle in Earth system (rad). Not used in this method.
 *
 * Returns the energy dispersion in units of s^-1 MeV^-1 for a given observed energy,
 * true photon energy and offset angle.
 ***************************************************************************/
double GCTAEdisp2D::operator()(const double& logEobs, 
                              const double& logEsrc, 
                              const double& theta, 
                              const double& phi,
                              const double& zenith,
                              const double& azimuth) const
{
    // Initalize edisp
    double edisp = 0.0;
    
    // Compute edisp
    edisp = m_edisp(0, logEsrc, Eobs / Esrc, theta);
    
    // Return
    return edisp;
}

/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear instance
 *
 * This method properly resets the object to an initial state.
 ***************************************************************************/
void GCTAEdisp2D::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GCTAEdisp::free_members();

    // Initialise members
    this->GCTAEdisp::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
 *
 * @return Deep copy of energy dispersion instance.
 ***************************************************************************/
GCTAEdisp2D* GCTAEdisp2D::clone(void) const
{
    return new GCTAEdisp2D(*this);
}


/***********************************************************************//**
 * @brief Load energy dispersion from FITS file
 *
 * @param[in] filename FITS file.
 *
 * This method loads the energy dispersion information from a FITS file.
 ***************************************************************************/
void GCTAEdisp2D::load(const std::string& filename)
{
    // Open FITS file
    GFits fits(filename);

    // Read energy dispersion from file
    read(fits);

    // Close FITS file
    fits.close();

    // Store filename
    m_filename = filename;

    // Return
    return;
}

/***********************************************************************//**
 * @brief Return response table
 *
 * @return String containing the class name ("GCTAEdisp2D").
 ***************************************************************************/
const GCTAResponseTable& GCTAEdisp2D::table(void) const
{
    return m_edisp;
}

/***********************************************************************//**
 * @brief Set response table
 *
 * @param[in] Response table.
 ***************************************************************************/
void GCTAEdisp2D::table(const GCTAResponseTable& table)
{
    m_edisp = table;
}

/***********************************************************************//**
 * @brief Write CTA energy dispersion table into FITS binary table object.
 *
 * @param[in] hdu FITS binary table.
 *
 * @todo Add necessary keywords.
 ***************************************************************************/
void GCTAEdisp2D::write(GFitsBinTable& hdu) const
{
    // Write background table
    m_edisp.write(hdu);

    // Return
    return;
}

/***********************************************************************//**
 * @brief Save energy dispersion table into FITS file
 *
 * @param[in] filename Energy dispersion table FITS file name.
 * @param[in] clobber Overwrite existing file? (true=yes)
 *
 * Save the energy dispersion table into a FITS file.
 ***************************************************************************/
void GCTAEdisp2D::save(const std::string& filename, const bool& clobber) const
{
    // Create binary table
    GFitsBinTable table;
    table.extname("ENERGY DISPERSION");

    // Write the Energy dispersion table
    write(table);

    // Create FITS file, append table, and write into the file
    GFits fits;
    fits.append(table);
    fits.saveto(filename, clobber);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read energy dispersion from FITS file
 *
 * @param[in] fits FITS file pointer.
 *
 * Reads the energy dispersion form the FITS file extension "ENERGY DISPERSION".
 * The data are stored in m_edisp which is of type GCTAResponseTable. The
 * energy axis will be set to log10, the offset angle axis to radians.
 ***************************************************************************/
void GCTAEdisp2D::read(const GFits& fits)
{
    // Clear response table
    m_edisp.clear();

    // Get energy dispersion table
    const GFitsTable& table = *fits.table("ENERGY DISPERSION");

    // Read energy dispersion table
    m_edisp.read(table);

    int etrue_size = m_edisp.axis(0);
    int migra_size = m_edisp.axis(1);
    int theta_size = m_edisp.axis(2);
    
    // Normalize vectors of migration matrix for all etrue and theta
    for (int i_etrue = 0; i_etrue < etrue_size; ++i_etrue) {
        double EobsOnEtrue  = 0.5 * (m_edisp.axis_hi(1, i_etrue) + m_edisp.axis_lo(1, i_etrue));
        for (int i_theta = 0; i_theta < theta_size; ++i_theta) {
            // Compute sum
            double sum = 0.0;
            for (int i_migra = 0; i_migra < migra_size; ++i_migra) {
                // Compute delta(Eobs/Etrue)
                double delta = m_edisp.axis_hi(1, i_migra) - m_edisp.axis_lo(1, i_migra);
                sum += m_edisp(0, i_etrue + (i_migra + i_theta*migra_size)*etrue_size) * delta / EobsOnEtrue;
            }
            if (sum > 0.0) {
                for (int i_migra = 0; i_migra < migra_size; ++i_migra) {
                    // Normalize by sum
                    m_edisp(0, i_etrue + (i_migra + i_theta*migra_size)*etrue_size) /= sum;
                }
            }
        }
    }
    
    // Set true energy axis to logarithmic scale
    m_edisp.axis_log10(0);
    
    // Set migration energy ratio axis to logarithmic scale
    m_edisp.axis_log10(1);

    // Set offset angle axis to radians
    m_edisp.axis_radians(2);

    // Return
    return;
}

/***********************************************************************//**
 * @brief Simulate energy dispersion
 *
 * @param[in] ran Random number generator.
 * @param[in] logE Log10 of the true photon energy (TeV).
 * @param[in] theta Offset angle in camera system (rad). Not used.
 * @param[in] phi Azimuth angle in camera system (rad). Not used.
 * @param[in] zenith Zenith angle in Earth system (rad). Not used.
 * @param[in] azimuth Azimuth angle in Earth system (rad). Not used.
 *
 * Draws observed energy value from a normal distribution of width
 * m_par_sigma around @p logE.
 ***************************************************************************/
GEnergy GCTAEdisp2D::mc(GRan&         ran,
                        const double& logEsrc,
                        const double& theta,
                        const double& phi,
                        const double& zenith,
                        const double& azimuth) const
{
    // Vector representing cumulative probability distribution
    std::vector<std::pair<double, double> > cumul;

    const double deltaLogE = 0.01;
    double logEobs = 0.01;
    double sum = 0.0;
    // Loop through Eobs
    for(int i = 0; i < 3.0 / deltaLogE; ++i) {
        // Increment Eobs
        logEobs += deltaLogE;
        
        // Compute cumulated probability
        sum += GCTAEdisp2D::operator()(logEobs, logEsrc, theta) * deltaLogE;
        
        // Create pair containing Eobs/Esrc and cumulated probability
        std::pair<double, double> pair(std::exp((logEobs - logEsrc) * std::log(10.0)), sum);
        
        // Add to vector
        cumul.push_back(pair);
    }
    
    // Draw random number between 0 and 1
    double p = ran.uniform();
    
    // Find right index
    int index = 0;
    while(index < 3.0 / deltaLogE - 1.0 && cumul[index+1].second < p) {
        index++;
    } // index found
    
    // Compute result Eobs/Esrc value
    double result =   (cumul[index+1].second - p)*cumul[index].first
                    + (p - cumul[index].second)*cumul[index+1].first;
    result       *=   (cumul[index+1].second - cumul[index].second);
    
    
    // Final result
    double logEnergy = std::log10(result) + logEsrc;
    
    // Set energy
    GEnergy energy;
    energy.log10TeV(logEnergy);

    // Return energy
    return energy;
}


/***********************************************************************//**
 * @brief Return observed energy interval that contains the energy dispersion.
 *
 * @param[in] logEsrc Log10 of the true photon energy (TeV).
 * @param[in] theta Offset angle in camera system (rad). Not used.
 * @param[in] phi Azimuth angle in camera system (rad). Not used.
 * @param[in] zenith Zenith angle in Earth system (rad). Not used.
 * @param[in] azimuth Azimuth angle in Earth system (rad). Not used.
 *
 * Returns the band of observed energies outside of which the energy
 * dispersion becomes negligible for a given true energy @p logEsrc. This
 * band is set to \f$\pm 5 \times \sigma\f$, where \f$\sigma\f$ is the
 * Gaussian width of the energy dispersion.
 ***************************************************************************/
GEbounds GCTAEdisp2D::ebounds_obs(const double& logEsrc,
                                  const double& theta,
                                  const double& phi,
                                  const double& zenith,
                                  const double& azimuth) const
{
    // Set energy band constant
    const double number_of_sigmas = 5.0;

    // Update the parameter cache
    update(m_par_logEobs, logEsrc, theta);

    // Compute energy boundaries
    GEnergy emin;
    GEnergy emax;
    emin.log10TeV(logEsrc - number_of_sigmas * m_eres);
    emax.log10TeV(logEsrc + number_of_sigmas * m_eres);

    // Return energy boundaries
    return (GEbounds(emin, emax));
}


/***********************************************************************//**
 * @brief Return true energy interval that contains the energy dispersion.
 *
 * @param[in] logEobs Log10 of the observed event energy (TeV).
 * @param[in] theta Offset angle in camera system (rad). Not used.
 * @param[in] phi Azimuth angle in camera system (rad). Not used.
 * @param[in] zenith Zenith angle in Earth system (rad). Not used.
 * @param[in] azimuth Azimuth angle in Earth system (rad). Not used.
 *
 * Returns the band of true photon energies outside of which the energy
 * dispersion becomes negligible for a given observed energy @p logEobs. This
 * band is set to \f$\pm 5 \times \sigma\f$, where \f$\sigma\f$ is the
 * Gaussian width of the energy dispersion.
 ***************************************************************************/
GEbounds GCTAEdisp2D::ebounds_src(const double& logEobs,
                                  const double& theta,
                                  const double& phi,
                                  const double& zenith,
                                  const double& azimuth) const
{
    // Set energy band constant
    const double number_of_sigmas = 5.0;

    // Update the parameter cache
    update(logEobs, m_par_logEsrc, theta);

    // Compute energy boundaries
    GEnergy emin;
    GEnergy emax;
    emin.log10TeV(logEobs - number_of_sigmas * m_eres);
    emax.log10TeV(logEobs + number_of_sigmas * m_eres);

    // Return energy boundaries
    return (GEbounds(emin, emax));
}

/***********************************************************************//**
 * @brief Print energy dispersion information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing energy dispersion information.
 ***************************************************************************/
std::string GCTAEdisp2D::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Compute energy boundaries in TeV
        double emin = m_edisp.axis_lo(0,0);
        double emax = m_edisp.axis_hi(0,m_edisp.axis(0)-1);

        // Compute offset angle boundaries in deg
        double omin = m_edisp.axis_lo(1,0);
        double omax = m_edisp.axis_hi(1,m_edisp.axis(1)-1);

        // Append header
        result.append("=== GCTAEdisp2D ===");

        // Append information
        result.append("\n"+gammalib::parformat("Filename")+m_filename);
        result.append("\n"+gammalib::parformat("Number of energy bins") +
                      gammalib::str(m_edisp.axis(0)));
        result.append("\n"+gammalib::parformat("Number of offset bins") +
                      gammalib::str(m_edisp.axis(1)));
        result.append("\n"+gammalib::parformat("Log10(Energy) range"));
        result.append(gammalib::str(emin)+" - "+gammalib::str(emax)+" TeV");
        result.append("\n"+gammalib::parformat("Offset angle range"));
        result.append(gammalib::str(omin)+" - "+gammalib::str(omax)+" deg");

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
void GCTAEdisp2D::init_members(void)
{
    // Initialise members
    m_filename.clear();
    m_edisp.clear();
    m_par_logEobs = -1.0e30;
    m_par_logEsrc = -1.0e30;
    m_par_theta = -1.0;
    m_norm = 1.0;
    m_EobsOnEsrc = 0.0;
    m_eres = 0.0;
    m_ebias = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] edisp Energy dispersion.
 ***************************************************************************/
void GCTAEdisp2D::copy_members(const GCTAEdisp2D& edisp)
{
    // Copy members
    m_filename = edisp.m_filename;
    m_edisp    = edisp.m_edisp;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAEdisp2D::free_members(void)
{
    // Return
    return;
}

/***********************************************************************//**
 * @brief Update Energy dispersion parameter cache
 *
 * @param[in] logEobs Log10 of the observed photon energy (TeV).
 * @param[in] logEsrc Log10 of the true photon energy (TeV).
 * @param[in] theta Offset angle in camera system (rad).
 *
 * This method updates the Energy dispersion parameter cache.
 ***************************************************************************/
void GCTAEdisp2D::update(const double& logEobs, const double& logEsrc, const double& theta) const
{
    /*
    // Only compute Edisp parameters if arguments have changed
    if (logEobs != m_par_logEobs || logEsrc != m_par_logEsrc || theta != m_par_theta) {

        // Save parameters
        m_par_logEobs  = logEobs;
        m_par_logEsrc  = logEsrc;
        m_par_theta = theta;

        // Compute linear Eobs/Esrc ratio
        m_EobsOnEsrc = std::exp( (logEobs - logEsrc) * std::log(10));
        
        
        
        
        // Get ERES and EBIAS values from response table
        m_eres = m_edisp(0, logEsrc, theta);
        m_ebias = m_edisp(1, logEsrc, theta);

        // Compute global normalization parameter
        m_norm = (m_eres > 0.0) ? gammalib::inv_sqrt2pi / m_eres : 0.0;

    }
    */
    // Return
    return;
}
