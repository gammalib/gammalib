/***************************************************************************
 *            GCTAMeanPsf.hpp - CTA mean point spread function class       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Chia-Chun Lu                                     *
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
 * @file GCTAMeanPsf.hpp
 * @brief CTA mean point spread function class definition
 * @author Chia-Chun Lu
 */

#ifndef GCTAMEANPSF_HPP
#define GCTAMEANPSF_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GFits.hpp"
#include "GSkymap.hpp"
#include "GObservations.hpp"
#include "GNodeArray.hpp"

/***********************************************************************//**
 * @class GCTAMeanPsf
 *
 * @brief class for the CTA point spread function
 *
 * This class implements class for the CTA point spread
 * function.
 ***************************************************************************/
class GCTAMeanPsf : public GBase {

public:
   
    // Constructors and destructors
    GCTAMeanPsf(void);
    GCTAMeanPsf(const GCTAMeanPsf& psf);
    GCTAMeanPsf(const GObservations& obs, 
		const double& x, const double& y, 
		const double& dx, const double& dy,
		const int& nx, const int& ny,
		const double& emin, const double& emax, const int& nebins,
		const double& min, const double& max, const int& nbins);
    virtual ~GCTAMeanPsf(void);

    // Operators
    GCTAMeanPsf& operator=(const GCTAMeanPsf& psf);

    // Methods
    void         clear(void);
    GCTAMeanPsf* clone(void) const;
    void         load(const std::string& filename);
    void         write(GFits& file) const;
    void         save(const std::string& filename, 
		      const bool& clobber) const;
    std::string  print(const GChatter& chatter = NORMAL) const;
    const GSkymap& map(void) const;
    const GEbounds& ebounds(void) const;
    const GNodeArray& deltas(void) const;

protected:
    // Methods
    void init_members(void);
    void copy_members(const GCTAMeanPsf& psf);
    void free_members(void);
    void set_psfcube(void);
    // Data
    GObservations m_obs; //!< Observation container
    int m_nbins; //!< number of delta bins
    int m_nebins; //!< number of energy bins
    GSkymap m_cube; //!< Average PSF cube
    GEbounds m_ebounds;  //!< Energy bounds for the PSF cube
    GNodeArray m_deltas; //!< delta bins for the PSF cube
   
};
/***********************************************************************//**
 * @brief Return psf cube sky map
 *
 * @return psf cube sky map.
 *
 * The GCTAMeanPsf represents the psf cube as a sky map. This methods
 * returns the sky map that is stored internally by GCTAMeanPsf as psf
 * cube.
 ***************************************************************************/
inline
const GSkymap& GCTAMeanPsf::map(void) const
{
    return (m_cube);
}

/***********************************************************************//**
 * @brief Return energy boundaries
 *
 * @return Energy boundaris
 *
 ***************************************************************************/
inline
const GEbounds& GCTAMeanPsf::ebounds(void) const
{
    return (m_ebounds);
}

/***********************************************************************//**
 * @brief Return deltas nodes
 *
 * @return deltas
 *
 ***************************************************************************/
inline
const GNodeArray& GCTAMeanPsf::deltas(void) const
{
    return (m_deltas);
}

#endif /* GCTAMEANPSF_HPP */
