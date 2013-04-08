/***************************************************************************
 *             GCTAEventCube.hpp - CTA event bin container class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2013 by Juergen Knoedlseder                         *
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
 * @file GCTAEventCube.hpp
 * @brief CTA event bin container class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTAEVENTCUBE_HPP
#define GCTAEVENTCUBE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GEventCube.hpp"
#include "GCTAEventBin.hpp"
#include "GSkymap.hpp"
#include "GEbounds.hpp"
#include "GGti.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GCTAInstDir.hpp"
#include "GFitsTable.hpp"
#include "GFitsImage.hpp"


/***********************************************************************//**
 * @class GCTAEventCube
 *
 * @brief CTA event bin container class
 *
 * This class is a container class for CTA event bins.
 ***************************************************************************/
class GCTAEventCube : public GEventCube {

public:
    // Constructors and destructors
    GCTAEventCube(void);
    explicit GCTAEventCube(const GSkymap& map, const GEbounds& ebds, const GGti& gti);
    GCTAEventCube(const GCTAEventCube& cube);
    virtual ~GCTAEventCube(void);

    // Operators
    virtual GCTAEventCube&      operator=(const GCTAEventCube& cube);
    virtual GCTAEventBin*       operator[](const int& index);
    virtual const GCTAEventBin* operator[](const int& index) const;

    // Implemented pure virtual base class methods
    virtual void           clear(void);
    virtual GCTAEventCube* clone(void) const;
    virtual int            size(void) const;
    virtual int            dim(void) const;
    virtual int            naxis(int axis) const;
    virtual void           load(const std::string& filename);
    virtual void           save(const std::string& filename, bool clobber = false) const;
    virtual void           read(const GFits& file);
    virtual void           write(GFits& file) const;
    virtual int            number(void) const;
    virtual std::string    print(const GChatter& chatter = NORMAL) const;

    // Other methods
    void                   map(const GSkymap& map);
    const GSkymap&         map(void) const { return m_map; }
    int                    nx(void) const { return m_map.nx(); }
    int                    ny(void) const { return m_map.ny(); }
    int                    npix(void) const { return m_map.npix(); }
    int                    ebins(void) const { return m_map.nmaps(); }

protected:
    // Protected methods
    void         init_members(void);
    void         copy_members(const GCTAEventCube& cube);
    void         free_members(void);
    void         read_cntmap(const GFitsImage* hdu);
    void         read_ebds(const GFitsTable* hdu);
    void         read_gti(const GFitsTable* hdu);
    void         set_directions(void);
    virtual void set_energies(void);
    virtual void set_times(void);
    void         set_bin(const int& index);

    // Protected members
    GSkymap                  m_map;        //!< Counts map stored as sky map
    GCTAEventBin             m_bin;        //!< Actual event bin
    GTime                    m_time;       //!< Event cube mean time
    std::vector<GCTAInstDir> m_dirs;       //!< Array of event directions
    std::vector<double>      m_omega;      //!< Array of solid angles (sr)
    std::vector<GEnergy>     m_energies;   //!< Array of log mean energies
    std::vector<GEnergy>     m_ewidth;     //!< Array of energy bin widths
    double                   m_ontime;     //!< Event cube ontime (sec)
};

#endif /* GCTAEVENTCUBE_HPP */
