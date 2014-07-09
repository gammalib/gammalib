/***************************************************************************
 *            GCTAExposure.hpp - CTA mean point spread function class       *
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
 * @file GCTAExposure.hpp
 * @brief CTA mean point spread function class definition
 * @author Chia-Chun Lu
 */

#ifndef GCTAEXPOSURE_HPP
#define GCTAEXPOSURE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GFits.hpp"
#include "GSkymap.hpp"
#include "GObservations.hpp"
#include "GNodeArray.hpp"

/***********************************************************************//**
 * @class GCTAExposure
 *
 * @brief Class of CTA exposure cube
 *
 * This class implements the CTA exposure cube
 *
 ***************************************************************************/
class GCTAExposure : public GBase {

public:
   
    // Constructors and destructors
    GCTAExposure(void);
    GCTAExposure(const GCTAExposure& exp);
    GCTAExposure(const GObservations& obs, 
		const double& x, const double& y, 
		const double& dx, const double& dy,
		const int& nx, const int& ny,
		const double& emin, const double& emax, const int& nebins);
    virtual ~GCTAExposure(void);

    // Operators
    GCTAExposure& operator=(const GCTAExposure& exp);

    // Methods
    void         clear(void);
    GCTAExposure* clone(void) const;
    void         load(const std::string& filename);
    void         write(GFits& file) const;
    void         save(const std::string& filename, 
		      const bool& clobber) const;
    std::string  print(const GChatter& chatter = NORMAL) const;
    const GSkymap& map(void) const;
    const GEbounds& ebounds(void) const;

protected:
    // Methods
    void init_members(void);
    void copy_members(const GCTAExposure& exp);
    void free_members(void);
    void set_expcube(void);

    // Data
    GObservations m_obs; //!< Observation container
    int m_nbins; //!< number of delta bins
    int m_nebins; //!< number of energy bins
    GSkymap m_cube; //!< Average Exposure cube
    GEbounds m_ebounds;  //!< Energy bounds for the Exposure cube
   
};
/***********************************************************************//**
 * @brief Return Exposure cube sky map
 *
 * @return Exposure cube sky map.
 *
 * The GCTAExposure represents the Exposure cube as a sky map. This methods
 * returns the sky map that is stored internally by GCTAExposure as Exposure
 * cube.
 ***************************************************************************/
inline
const GSkymap& GCTAExposure::map(void) const
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
const GEbounds& GCTAExposure::ebounds(void) const
{
    return (m_ebounds);
}



#endif /* GCTAEXPOSURE_HPP */
