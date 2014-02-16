/***************************************************************************
 *            GCTAEdispRmf.cpp - CTA RMF energy dispersion class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Christoph Deil & Ellis Owen                      *
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
 * @file GCTAEdispRmf.cpp
 * @brief CTA RMF energy dispersion class implementation
 * @author Christoph Deil & Ellis Owen
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include "GTools.hpp"
#include "GMath.hpp"
#include "GEnergy.hpp"
#include "GCTAEdispRmf.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_LOAD                             "GCTAEdispRmf::load(std::string&)"

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
GCTAEdispRmf::GCTAEdispRmf(void) : GCTAEdisp()
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief File constructor
 *
 * @param[in] filename Redistribution Matrix File name.
 *
 * Construct instance by loading the Redistribution Matrix File.
 ***************************************************************************/
GCTAEdispRmf::GCTAEdispRmf(const std::string& filename) : GCTAEdisp()
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
 * @param[in] edisp Energy dispersion
 ***************************************************************************/
GCTAEdispRmf::GCTAEdispRmf(const GCTAEdispRmf& edisp) : GCTAEdisp(edisp)
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
GCTAEdispRmf::~GCTAEdispRmf(void)
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
 * @param[in] edisp Energy dispersion
 * @return Energy dispersion
 ***************************************************************************/
GCTAEdispRmf& GCTAEdispRmf::operator=(const GCTAEdispRmf& edisp)
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
 * @brief Return energy dispersion.
 *
 * @param[in] logEobs log10 of the observed photon energy (TeV).
 * @param[in] logEsrc log10 of the true photon energy (TeV).
 * @param[in] theta Offset angle in camera system (rad). Not used.
 * @param[in] phi Azimuth angle in camera system (rad). Not used.
 * @param[in] zenith Zenith angle in Earth system (rad). Not used.
 * @param[in] azimuth Azimuth angle in Earth system (rad). Not used.
 *
 * Returns the energy resolution, i.e. the probability density in observed
 * photon energy at a given (log10(E_src), log10(E_obs)).
 * To be precise: energy dispersion = dP / d(log10(E_obs)).
 * 
 * @todo Implement interpolation method
 * @todo So far the operator returns a matrix element and not a density!!!
 ***************************************************************************/
double GCTAEdispRmf::operator()(const double& logEobs,
                                const double& logEsrc,
                                const double& theta,
                                const double& phi,
                                const double& zenith,
                                const double& azimuth) const
{
    // Set true energy
    GEnergy etrue;
    etrue.log10TeV(logEsrc);

    // Set measured energy
    GEnergy emeasured;
    emeasured.log10TeV(logEobs);

    // Get indices for observed and true energy 
	int itrue     = m_rmf.etrue().index(etrue);
	int imeasured = m_rmf.emeasured().index(emeasured);

    // Extract matrix element
	double edisp = m_rmf.matrix()(itrue, imeasured);
    
    // Return energy dispersion
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
void GCTAEdispRmf::clear(void)
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
 * @return Deep copy of instance.
 ***************************************************************************/
GCTAEdispRmf* GCTAEdispRmf::clone(void) const
{
    return new GCTAEdispRmf(*this);
}


/***********************************************************************//**
 * @brief Load energy dispersion from RMF file
 *
 * @param[in] filename of RMF file.
 *
 * This method loads the energy dispersion information from an RMF file.
 ***************************************************************************/
void GCTAEdispRmf::load(const std::string& filename)
{
    // Load RMF file
    m_rmf.load(filename);

    // Store the filename
    m_filename = filename;

    // Set Monte Carlo cache
    set_mc_cache();

	// Return
    return;
}


/***********************************************************************//**
 * @brief Simulate energy dispersion
 *
 * @param[in] ran Random number generator.
 * @param[in] logEsrc Log10 of the true photon energy (TeV).
 * @param[in] theta Offset angle in camera system (rad). Not used.
 * @param[in] phi Azimuth angle in camera system (rad). Not used.
 * @param[in] zenith Zenith angle in Earth system (rad). Not used.
 * @param[in] azimuth Azimuth angle in Earth system (rad). Not used.
 *
 * Draws observed energy value from RMF matrix.
 ***************************************************************************/
GEnergy GCTAEdispRmf::mc(GRan&         ran,
                         const double& logEsrc,
                         const double& theta,
                         const double& phi,
                         const double& zenith,
                         const double& azimuth) const
{
    // Set true energy
    GEnergy energy;
    energy.log10TeV(logEsrc);

    // Determine true energy index
    int itrue = m_rmf.etrue().index(energy);

    // Continue only if the true energy index lies within a true energy
    // boundary of the redistribution matrix file
    if (itrue != -1) {

        // Get offset in measured energy. Continue only if offset is valid
        int offset = m_mc_measured_start[itrue];
        if (offset != -1) {

            // Determine measured energy index from Monte-Carlo cache
            int imeasured = ran.cdf(m_mc_measured_cdf[itrue]) + offset;

            // Get log10 energy minimum and bin width
            double emin   = m_rmf.emeasured().emin(imeasured).log10TeV();
            double ewidth = m_rmf.emeasured().emax(imeasured).log10TeV() - emin;

            // Draw a random energy from this interval
            double e = emin + ewidth * ran.uniform();

            // Set interval
            energy.log10TeV(e);
        
        } // endif: offset was valid

    } // endif: there was redistribution information

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
 * dispersion becomes negligible for a given true energy @p logEsrc.
 ***************************************************************************/
GEbounds GCTAEdispRmf::ebounds_obs(const double& logEsrc,
                                   const double& theta,
                                   const double& phi,
                                   const double& zenith,
                                   const double& azimuth) const
{
    // Set true energy
    GEnergy etrue;
    etrue.log10TeV(logEsrc);

    // Return measured energy boundaries
    return (m_rmf.emeasured(etrue));
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
 * dispersion becomes negligible for a given observed energy @p logEobs.
 ***************************************************************************/
GEbounds GCTAEdispRmf::ebounds_src(const double& logEobs,
                                   const double& theta,
                                   const double& phi,
                                   const double& zenith,
                                   const double& azimuth) const
{
    // Set measured energy
    GEnergy emeasured;
    emeasured.log10TeV(logEobs);

    // Return true energy boundaries
    return (m_rmf.etrue(emeasured));
}


/***********************************************************************//**
 * @brief Print RMF information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing energy dispersion information.
 ***************************************************************************/
std::string GCTAEdispRmf::print(const GChatter& chatter) const
{ 
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCTAEdispRmf ===");

        // Append energy boundary information
        result.append("\n"+gammalib::parformat("Filename")+m_filename);
        result.append("\n"+gammalib::parformat("Number of true energy bins"));
        result.append(gammalib::str(m_rmf.etrue().size()));
        result.append("\n"+gammalib::parformat("Number of measured bins"));
        result.append(gammalib::str(m_rmf.emeasured().size()));
        result.append("\n"+gammalib::parformat("True energy range"));
        result.append(m_rmf.etrue().emin().print());
        result.append(" - ");
        result.append(m_rmf.etrue().emax().print());
        result.append("\n"+gammalib::parformat("Measured energy range"));
        result.append(m_rmf.emeasured().emin().print());
        result.append(" - ");
        result.append(m_rmf.emeasured().emax().print());

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
void GCTAEdispRmf::init_members(void)
{
    // Initialise members
    m_filename.clear();
    m_rmf.clear();

    // Initialise cache
    m_mc_measured_start.clear();
    m_mc_measured_cdf.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] edisp Energy dispersion
 ***************************************************************************/
void GCTAEdispRmf::copy_members(const GCTAEdispRmf& edisp)
{
    // Copy members
    m_filename = edisp.m_filename;
    m_rmf      = edisp.m_rmf;

    // Copy cache
    m_mc_measured_start = edisp.m_mc_measured_start;
    m_mc_measured_cdf   = edisp.m_mc_measured_cdf;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAEdispRmf::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Update the Monte-Carlo cache values
 *
 * Sets the following private class members:
 *
 *      m_mc_measured_start: start index in emeasured boundaries for each true energy
 *      m_mc_measured_stop:  stop index in emeasured boundaries for each true energy
 *      m_mc_measured_cdf:   CDF for each true energy
 *
 ***************************************************************************/
void GCTAEdispRmf::set_mc_cache(void)
{
    // Clear MC cache
    m_mc_measured_start.clear();
    m_mc_measured_cdf.clear();

    // Reserve some space
    m_mc_measured_start.reserve(m_rmf.ntrue());
    m_mc_measured_cdf.reserve(m_rmf.ntrue());

    // Loop over true energies
    for (int itrue = 0; itrue < m_rmf.ntrue(); ++itrue) {

        // Determine first non-zero measured energy
        int imeasured_start = -1;
        for (int i = 0; i < m_rmf.nmeasured(); ++i) {
            if (m_rmf(itrue, i) > 0.0) {
                imeasured_start = i;
                break;
            }
        }
        m_mc_measured_start.push_back(imeasured_start);

        // Determine last non-zero measured energy
        int imeasured_stop = -1;
        for (int i = m_rmf.nmeasured()-1; i >= 0; --i) {
            if (m_rmf(itrue, i) > 0.0) {
                imeasured_stop = i;
                break;
            }
        }

        // Determine number of elements
        int num = (imeasured_stop - imeasured_start + 1);
        if (num < 0) {
            num = 0;
        }

        // Allocate vector
        GVector vector(num+1);

        // If there are measured energies for this true energy then build
        // now the CDF ...
        if (imeasured_start != -1 && imeasured_stop != -1) {
            double sum = 0.0;
            vector[0] = 0.0;
            for (int i = imeasured_start, k = 1; i <= imeasured_stop; ++i, ++k) {
                sum      += m_rmf(itrue, i);
                vector[k] = sum;
            }
            if (sum > 0.0) {
                for (int k = 0; k < num; ++k) {
                    vector[k] /= sum;
                }
            }
        }

        // Make sure that last pixel in the cache is >1
        vector[num] = 1.0001;

        // Push vector on cache
        m_mc_measured_cdf.push_back(vector);
    
    } // endfor: looped over all true energies

    // Return
	return;
}
