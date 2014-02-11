/***************************************************************************
 *  GCTAEdispRMF.cpp - CTA RMF energy dispersion class *
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
 * @file GCTAEdispRMF.cpp
 * @brief CTA RMF energy dispersion class implementation
 * @author Christoph Deil & Ellis Owen
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdio>             // std::fopen, std::fgets, and std::fclose
#include <cmath>
#include "GTools.hpp"
#include "GRmf.hpp"
#include "GMath.hpp"
#include "GIntegral.hpp"
#include "GFitsTable.hpp"
#include "GFitsTableCol.hpp"
#include "GCTAEdispRMF.hpp"
#include "GCTAResponse.hpp"
#include "GCTAResponse_helpers.hpp"
#include "GCTAException.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_LOAD                       "GCTAEdispRMF::load(std::string&)"

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
GCTAEdispRMF::GCTAEdispRMF(void) : GCTAEdisp()
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief File constructor
 *
 * @param[in] filename Performance table file name.
 *
 * Construct instance by loading the energy dispersion information from
 * an ASCII performance table.
 ***************************************************************************/
GCTAEdispRMF::GCTAEdispRMF(const std::string& filename) :
                    GCTAEdisp()
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
GCTAEdispRMF::GCTAEdispRMF(const GCTAEdispRMF& edisp) :
                    GCTAEdisp(edisp)
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
GCTAEdispRMF::~GCTAEdispRMF(void)
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
GCTAEdispRMF& GCTAEdispRMF::operator=(const GCTAEdispRMF& edisp)
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
 * Evaluates
 *
 * \f[
 * S(E) = \frac{1}{\sqrt{2\pi}m\_sigma}
 *        \exp(\frac{-(logEobs-logEsrc)^2}{2 m\_sigma^2})
 * \f]
 ***************************************************************************/
double GCTAEdispRMF::operator()(const double& logEobs,
                                      const double& logEsrc,
                                      const double& theta,
                                      const double& phi,
                                      const double& zenith,
                                      const double& azimuth) const
{

	int i_src = m_ebounds_src.index(GEnergy(std::pow(10, logEsrc), "TeV"));
	int i_obs = m_ebounds_src.index(GEnergy(std::pow(10, logEobs), "TeV"));
	double edisp = m_matrix(i_src, i_obs);
    
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
void GCTAEdispRMF::clear(void)
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
GCTAEdispRMF* GCTAEdispRMF::clone(void) const
{
    return new GCTAEdispRMF(*this);
}


/***********************************************************************//**
 * @brief Load energy dispersion from performance table
 *
 * @param[in] filename Performance table file name.
 *
 * @exception GCTAExceptionHandler::file_open_error
 *            File could not be opened for read access.
 *
 * This method loads the energy dispersion information from an ASCII
 * performance table. The energy resolution is stored in the 5th column
 * of the performance table as RMS(ln(Eest/Etrue)). The method converts
 * this internally to a sigma value by multiplying the stored values by
 * 0.434294481903.
 ***************************************************************************/
void GCTAEdispRMF::load(const std::string& filename)
{
	GRmf rmf;
	rmf.filename();

	m_ebounds_src = rmf.etrue();
	m_ebounds_obs = rmf.emeasured();
	m_matrix = rmf.matrix();

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
GEnergy GCTAEdispRMF::mc(GRan&         ran,
                         const double& logE,
                         const double& theta,
                         const double& phi,
                         const double& zenith,
                         const double& azimuth) const
{
    // TODO
    // Return energy
    // return energy;
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
GEbounds GCTAEdispRMF::ebounds_obs(const double& logEsrc,
                                         const double& theta,
                                         const double& phi,
                                         const double& zenith,
                                         const double& azimuth) const
{
    // TODO
    // Return energy boundaries
    //return (GEbounds(emin, emax));
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
GEbounds GCTAEdispRMF::ebounds_src(const double& logEobs,
                                         const double& theta,
                                         const double& phi,
                                         const double& zenith,
                                         const double& azimuth) const
{
   // TODO

    // Return energy boundaries
    //return (GEbounds(emin, emax));
}


/***********************************************************************//**
 * @brief Print energy dispersion information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing energy dispersion information.
 ***************************************************************************/
std::string GCTAEdispRMF::print(const GChatter& chatter) const
{
	// TODO

    // Return result
    //return result;
}


/*==========================================================================
 =                                                                         =
 =                            Private methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GCTAEdispRMF::init_members(void)
{
    // Initialise members
    m_filename.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] edisp Energy dispersion
 ***************************************************************************/
void GCTAEdispRMF::copy_members(const GCTAEdispRMF& edisp)
{
    // Copy members
    m_filename  = edisp.m_filename;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAEdispRMF::free_members(void)
{
    // Return
    return;
}

