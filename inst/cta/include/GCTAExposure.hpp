/***************************************************************************
 *                GCTAExposure.hpp - CTA exposure cube class               *
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
 * @brief CTA exposure cube class definition
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
#include "GCTAObservation.hpp"


/***********************************************************************//**
 * @class GCTAExposure
 *
 * @brief CTA exposure cube class
 *
 * This class implements a CTA exposure cube which provides the average
 * exposure for binned analysis as function of sky position and log10
 * energy.
 ***************************************************************************/
class GCTAExposure : public GBase {

public:
   
    // Constructors and destructors
    GCTAExposure(void);
    GCTAExposure(const GCTAExposure& cube);
    GCTAExposure(const std::string&   wcs,
                 const std::string&   coords,
                 const double&        x,
                 const double&        y,
                 const double&        dx,
                 const double&        dy,
                 const int&           nx,
                 const int&           ny,
                 const GEbounds&      ebounds);
    virtual ~GCTAExposure(void);

    // Operators
    GCTAExposure& operator=(const GCTAExposure& exp);

    // Methods
    void            clear(void);
    GCTAExposure*   clone(void) const;
    void            set(const GCTAObservation& obs);
    void            fill(const GObservations& obs);
    const GSkymap&  cube(void) const;
    const GEbounds& ebounds(void) const;
    void            write(GFits& file) const;
    void            load(const std::string& filename);
    void            save(const std::string& filename,
                         const bool& clobber = false) const;
    std::string     print(const GChatter& chatter = NORMAL) const;

protected:
    // Methods
    void init_members(void);
    void copy_members(const GCTAExposure& exp);
    void free_members(void);
    void clear_cube(void);

    // Data
    GSkymap  m_cube;     //!< Average Exposure cube
    GEbounds m_ebounds;  //!< Energy bounds for the Exposure cube
};


/***********************************************************************//**
 * @brief Return exposure cube
 *
 * @return Exposure cube.
 *
 * Returns the GSkymap object that is used to store the exposure cube
 * information.
 ***************************************************************************/
inline
const GSkymap& GCTAExposure::cube(void) const
{
    return (m_cube);
}

/***********************************************************************//**
 * @brief Return energy boundaries
 *
 * @return Energy boundaries
 ***************************************************************************/
inline
const GEbounds& GCTAExposure::ebounds(void) const
{
    return (m_ebounds);
}

#endif /* GCTAEXPOSURE_HPP */
