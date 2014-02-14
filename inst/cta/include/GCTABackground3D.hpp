/***************************************************************************
 *              GCTABackground3D.hpp - CTA 3D background class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Juergen Knoedlseder                              *
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
 * @file GCTABackground3D.hpp
 * @brief CTA 3D background class definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTABACKGROUND3D_HPP
#define GCTABACKGROUND3D_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GCTABackground.hpp"
#include "GCTAResponseTable.hpp"

/* __ Forward declarations _______________________________________________ */
class GFits;


/***********************************************************************//**
 * @class GCTABackground3D
 *
 * @brief CTA 3D background class
 ***************************************************************************/
class GCTABackground3D : public GCTABackground {

public:
    // Constructors and destructors
    GCTABackground3D(void);
    explicit GCTABackground3D(const std::string& filename);
    GCTABackground3D(const GCTABackground3D& bgd);
    virtual ~GCTABackground3D(void);

    // Operators
    GCTABackground3D& operator=(const GCTABackground3D& bgd);

    // Implemented pure virtual methods
    void                       clear(void);
    GCTABackground3D*          clone(void) const;
    void                       load(const std::string& filename);
    std::string                filename(void) const;
    GCTAInstDir                mc(const GEnergy& energy,
                                  const GTime& time,
                                  GRan& ran) const;
    const GModelSpectralNodes& spectrum(void) const;
    std::string                print(const GChatter& chatter = NORMAL) const;

    // Methods
    void read(const GFits& file);
    
private:
    // Methods
    void init_members(void);
    void copy_members(const GCTABackground3D& bgd);
    void free_members(void);
    void init_mc_cache(void) const;

    // Members
    std::string         m_filename;    //!< Name of background response file
    GCTAResponseTable   m_background;  //!< Background response table

    // Monte Carlo cache
    mutable std::vector<double> m_mc_cache;    //!< Monte Carlo cache
    mutable GModelSpectralNodes m_mc_spectrum; //!< Response cube spectrum
};


/***********************************************************************//**
 * @brief Return filename
 *
 * @return Returns filename from which the background was loaded.
 ***************************************************************************/
inline
std::string GCTABackground3D::filename(void) const
{
    // Return filename
    return m_filename;
}


/***********************************************************************//**
 * @brief Get response cube spectrum
 *
 * @return Response cube spectrum.
 *
 * Returns the response cube spectrum.
 ***************************************************************************/
inline
const GModelSpectralNodes& GCTABackground3D::spectrum(void) const
{
    if (m_mc_spectrum.nodes() == 0) {
        init_mc_cache();
    }
    return (m_mc_spectrum);
}

#endif /* GCTABACKGROUND3D_HPP */
