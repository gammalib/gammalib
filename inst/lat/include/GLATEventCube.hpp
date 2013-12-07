/***************************************************************************
 *                 GLATEventCube.hpp - LAT event cube class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2013 by Juergen Knoedlseder                         *
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
 * @file GLATEventCube.hpp
 * @brief LAT event cube class interface definition.
 * @author Juergen Knoedlseder
 */

#ifndef GLATEVENTCUBE_HPP
#define GLATEVENTCUBE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GEventCube.hpp"
#include "GLATInstDir.hpp"
#include "GLATEventBin.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GFits.hpp"
#include "GFitsImage.hpp"
#include "GFitsTable.hpp"
#include "GSkymap.hpp"
#include "GNodeArray.hpp"


/***********************************************************************//**
 * @class GLATEventCube
 *
 * @brief LAT event cube class interface defintion
 ***************************************************************************/
class GLATEventCube : public GEventCube {

public:
    // Constructors and destructors
    GLATEventCube(void);
    GLATEventCube(const GLATEventCube& cube);
    virtual ~GLATEventCube(void);

    // Operators
    virtual GLATEventCube&      operator=(const GLATEventCube& cube);
    virtual GLATEventBin*       operator[](const int& index);
    virtual const GLATEventBin* operator[](const int& index) const;

    // Implemented pure virtual base class methods
    virtual void           clear(void);
    virtual GLATEventCube* clone(void) const;
    virtual int            size(void) const;
    virtual int            dim(void) const;
    virtual int            naxis(const int& axis) const;
    virtual void           load(const std::string& filename);
    virtual void           save(const std::string& filename,
                                const bool& clobber = false) const;
    virtual void           read(const GFits& file);
    virtual void           write(GFits& file) const;
    virtual int            number(void) const;
    virtual std::string    print(const GChatter& chatter = NORMAL) const;

    // Other methods
    void              time(const GTime& time) { m_time=time; }
    void              map(const GSkymap& map);
    void              enodes(const GNodeArray& enodes) { m_enodes=enodes; }
    void              ontime(const double& ontime) { m_ontime=ontime; }
    const GTime&      time(void) const { return m_time; }
    const GSkymap&    map(void) const { return m_map; }
    const GNodeArray& enodes(void) { return m_enodes; }
    const double&     ontime(void) const { return m_ontime; }
    int               nx(void) const { return m_map.nx(); }
    int               ny(void) const { return m_map.ny(); }
    int               npix(void) const { return m_map.npix(); }
    int               ebins(void) const { return m_map.nmaps(); }
    int               ndiffrsp(void) const { return m_srcmap.size(); }
    std::string       diffname(const int& index) const;
    GSkymap*          diffrsp(const int& index) const;
    double            maxrad(const GSkyDir& dir) const;

protected:
    // Protected methods
    void         init_members(void);
    void         copy_members(const GLATEventCube& cube);
    void         free_members(void);
    void         read_cntmap(const GFitsImage* hdu);
    void         read_srcmap(const GFitsImage* hdu);
    void         read_ebds(const GFitsTable* hdu);
    void         read_gti(const GFitsTable* hdu);
    void         set_directions(void);
    virtual void set_energies(void);
    virtual void set_times(void);
    void         set_bin(const int& index);

    // Protected data area
    GLATEventBin             m_bin;          //!< Actual energy bin
    GSkymap                  m_map;          //!< Counts map stored as sky map
    GTime                    m_time;         //!< Event cube mean time
    double                   m_ontime;       //!< Event cube ontime (sec)
    std::vector<GLATInstDir> m_dirs;         //!< Array of event directions
    std::vector<double>      m_omega;        //!< Array of solid angles (sr)
    std::vector<GEnergy>     m_energies;     //!< Array of log mean energies
    std::vector<GEnergy>     m_ewidth;       //!< Array of energy bin widths
    std::vector<GSkymap*>    m_srcmap;       //!< Pointers to source maps
    std::vector<std::string> m_srcmap_names; //!< Source map names
    GNodeArray               m_enodes;       //!< Energy nodes
};

#endif /* GLATEVENTCUBE_HPP */
