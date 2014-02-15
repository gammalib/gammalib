/***************************************************************************
 *            GCTAEdispRMF.cpp - CTA RMF energy dispersion class           *
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
#include "GRan.hpp"
#include "GEnergy.hpp"
#include "GMath.hpp"
#include "GVector.hpp"
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
GCTAEdispRMF::GCTAEdispRMF(const std::string& filename) : GCTAEdisp()
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
GCTAEdispRMF::GCTAEdispRMF(const GCTAEdispRMF& edisp) : GCTAEdisp(edisp)
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
 *
 * @todo Implement interpolation method
 ***************************************************************************/
double GCTAEdispRMF::operator()(const double& logEobs,
                                const double& logEsrc,
                                const double& theta,
                                const double& phi,
                                const double& zenith,
                                const double& azimuth) const
{
    // Get indices for observed and true energy 
	int i_obs = m_ebounds_obs.index(GEnergy(std::pow(10.0, logEobs), "TeV"));
	int i_src = m_ebounds_src.index(GEnergy(std::pow(10.0, logEsrc), "TeV"));

    // Extract matrix element
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
 * @brief Load energy dispersion from RMF file
 *
 * @param[in] filename of RMF file.
 *
 * @exception GCTAExceptionHandler::file_open_error
 *            File could not be opened for read access.
 *
 * This method loads the energy dispersion information from an RMF file.
 *
 * @todo Why not storing simply GRmf as a member of GCTAEdispRMF
 ***************************************************************************/
void GCTAEdispRMF::load(const std::string& filename)
{
    // Load RMF file
	GRmf rmf(filename);

    // Extract infomation
	m_ebounds_src = rmf.etrue();
	m_ebounds_obs = rmf.emeasured();
	m_matrix      = rmf.matrix();

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
 * Draws observed energy value from RMF matrix.
 *
 * @todo Are you sure this method works??? I think the vector index should
 * be extracted from logE. Also, never copy a vector, get references to it.
 * Another problem is that the method will always return the bin centres,
 * but we want to sample energies to infinite precision!!!
 ***************************************************************************/
GEnergy GCTAEdispRMF::mc(GRan&         ran,
                         const double& logE,
                         const double& theta,
                         const double& phi,
                         const double& zenith,
                         const double& azimuth) const
{
	// Random selection of a GVector from vector of GVectors
	double  u      = ran.uniform();
	GVector vector = m_cdf_cache[u]; // Does this work??? u=[0,1]

	// Random selection of element in GVector
	int ref = ran.cdf(vector);

   	// Convert index into an energy
    double logEobs = vector[ref];

    // Draw log observed energy in TeV
    GEnergy energy;
    energy.log10TeV(logEobs);

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
 *
 * @todo Convert logEsrc into index!!!
 ***************************************************************************/
GEbounds GCTAEdispRMF::ebounds_obs(const double& logEsrc,
                                   const double& theta,
                                   const double& phi,
                                   const double& zenith,
                                   const double& azimuth) const
{
    // WRONG!!! No random sampling needed here!!!
    GRan ran;
	int low = ran.cdf(m_mc_cache);
	GVector EVector = m_matrix.column(low);

    // Initialise energy boundaries
    GEbounds ebounds;

    // Determine vector column that corresponds to specified true energy
    GEnergy energy;
    energy.log10TeV(logEsrc);
    int column = m_ebounds_src.index(energy);
    if (column != -1) {
    
        // Determine first and last non-zero indices
        int i_emin = -1;
        int i_emax = -1;
        for (int row = 0; row < m_matrix.rows(); ++row) {
            if (m_matrix(row, column) > 0.0) {
                i_emin = row;
                break;
            }
        }
        for (int row = m_matrix.rows()-1; row >= 0; --row) {
            if (m_matrix(row, column) > 0.0) {
                i_emax = row;
                break;
            }
        }

        // Set energy boundaries if valid indices have been found
        if (i_emin != -1 && i_emax != -1) {
            ebounds = GEbounds(m_ebounds_obs.emin(i_emin),
                               m_ebounds_obs.emax(i_emax));
        }

    } // endif: information found for source energy

    // Return energy boundaries
    return ebounds;
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
GEbounds GCTAEdispRMF::ebounds_src(const double& logEobs,
                                   const double& theta,
                                   const double& phi,
                                   const double& zenith,
                                   const double& azimuth) const
{
    GRan ran;
	int low = ran.cdf(m_mc_cache);
	GVector EVector = m_matrix.row(low);

	int emin_index = EVector.first_nonzero();
	int emax_index = EVector.last_nonzero();

	double emin_val = EVector[emin_index];
	double emax_val = EVector[emax_index];

	GEnergy emin = GEnergy(emin_val, "MeV");
	GEnergy emax = GEnergy(emax_val, "MeV");

	// Return energy boundaries
	return (GEbounds(emin, emax));

	}

/***********************************************************************//**
 * @brief Calculates std::vector of cdf GVectors.
 *
 * Makes std::vector of cdf GVectors from RMF matrix columns and puts
 * result in m_cdf_cache in the cache.
 ***************************************************************************/
void GCTAEdispRMF::convert_cdf(void)
{

	GVector Row = m_matrix.row(0);
	int row_sz = Row.size();

	//Iterate over vectors
	for(int j = 0; j <= row_sz; ++j) {
		GVector Vector = m_matrix.column(j);
		//Get minimum and maximum non-zero indices for vector
		int min_index = Vector.first_nonzero();
		int max_index = Vector.last_nonzero();
		int size = Vector.size();
		//Construct new vector for these non-zero values
		GVector v1 = Vector.slice((Vector[0] + min_index), (Vector[size] - max_index));
		int v1_size = v1.size();
		//First element of the new vector
		double el_1 = v1[0];
		//Copy new vector for iteration
		GVector v_add = v1;


		//Compute the cdf from shortened vector
		double element = v_add[0];
		//Iterate over elements of each vector
		for (int i = 0; i <= size; ++i) {
			double element_i = v_add[i];
			double new_element_i = element + element_i;
			//Insert replacement element into the vector at position i
			v_add[i] = new_element_i;
			double element = element_i;
		}
		//Adds the GVector cdf to the std::vector of GVectors
		GVector vector_cdf = v_add;
		m_cdf_cache[j] = vector_cdf;
		}

	return;

}

/***********************************************************************//**
 * @brief Print RMF information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing energy dispersion information.
 *
 * TODO: Print something useful.
 ***************************************************************************/
std::string GCTAEdispRMF::print(const GChatter& chatter) const
{ // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

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

