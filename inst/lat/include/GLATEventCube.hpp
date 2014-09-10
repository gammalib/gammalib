/***************************************************************************
 *             GLATEventCube.hpp - Fermi/LAT event cube class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2014 by Juergen Knoedlseder                         *
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
 * @brief Fermi/LAT event cube class definition
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
 * @brief Fermi/LAT event cube class
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
    virtual std::string    classname(void) const;
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
    void              time(const GTime& time);
    void              map(const GSkymap& map);
    void              enodes(const GNodeArray& enodes);
    void              ontime(const double& ontime);
    const GTime&      time(void) const;
    const GSkymap&    map(void) const;
    const GNodeArray& enodes(void) const;
    const double&     ontime(void) const;
    int               nx(void) const;
    int               ny(void) const;
    int               npix(void) const;
    int               ebins(void) const;
    int               ndiffrsp(void) const;
    std::string       diffname(const int& index) const;
    GSkymap*          diffrsp(const int& index) const;
    double            maxrad(const GSkyDir& dir) const;

protected:
    // Protected methods
    void         init_members(void);
    void         copy_members(const GLATEventCube& cube);
    void         free_members(void);
    void         read_cntmap(const GFitsImage& hdu);
    void         read_srcmap(const GFitsImage& hdu);
    void         read_ebds(const GFitsTable& hdu);
    void         read_gti(const GFitsTable& hdu);
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
    std::vector<double>      m_solidangle;   //!< Array of solid angles (sr)
    std::vector<GEnergy>     m_energies;     //!< Array of log mean energies
    std::vector<GEnergy>     m_ewidth;       //!< Array of energy bin widths
    std::vector<GSkymap*>    m_srcmap;       //!< Pointers to source maps
    std::vector<std::string> m_srcmap_names; //!< Source map names
    GNodeArray               m_enodes;       //!< Energy nodes
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GLATEventCube").
 ***************************************************************************/
inline
std::string GLATEventCube::classname(void) const
{
    return ("GLATEventCube");
}


/***********************************************************************//**
 * @brief Set event cube mean time
 *
 * @param[in] time Event cube mean time.
 ***************************************************************************/
inline
void GLATEventCube::time(const GTime& time)
{
    m_time = time;
    return;
}


/***********************************************************************//**
 * @brief Set event cube energy nodes
 *
 * @param[in] enodes Energy nodes.
 ***************************************************************************/
inline
void GLATEventCube::enodes(const GNodeArray& enodes)
{
    m_enodes = enodes;
    return;
}


/***********************************************************************//**
 * @brief Set event cube ontime
 *
 * @param[in] ontime Ontime.
 ***************************************************************************/
inline
void GLATEventCube::ontime(const double& ontime)
{
    m_ontime = ontime;
    return;
}


/***********************************************************************//**
 * @brief Return event cube mean time
 *
 * @return Event cube mean time.
 ***************************************************************************/
inline
const GTime& GLATEventCube::time(void) const
{
    return m_time;
}


/***********************************************************************//**
 * @brief Return event cube sky map
 *
 * @return Sky map.
 ***************************************************************************/
inline
const GSkymap& GLATEventCube::map(void) const
{
    return m_map;
}


/***********************************************************************//**
 * @brief Return event cube energy nodes
 *
 * @return Energy nodes.
 ***************************************************************************/
inline
const GNodeArray& GLATEventCube::enodes(void) const
{
    return m_enodes;
}


/***********************************************************************//**
 * @brief Return event cube ontime
 *
 * @return Ontime.
 ***************************************************************************/
inline
const double& GLATEventCube::ontime(void) const
{
    return m_ontime;
}


/***********************************************************************//**
 * @brief Return number of bins in X direction
 *
 * @return Number of bins in X direction.
 ***************************************************************************/
inline
int GLATEventCube::nx(void) const
{
    return m_map.nx();
}


/***********************************************************************//**
 * @brief Return number of bins in Y direction
 *
 * @return Number of bins in Y direction.
 ***************************************************************************/
inline
int GLATEventCube::ny(void) const
{
    return m_map.ny();
}


/***********************************************************************//**
 * @brief Return number of pixels in event cube sky map
 *
 * @return Number of pixels in event cube sky map.
 ***************************************************************************/
inline
int GLATEventCube::npix(void) const
{
    return m_map.npix();
}


/***********************************************************************//**
 * @brief Return number of energy bins in event cube
 *
 * @return Number of energy bins in event cube.
 ***************************************************************************/
inline
int GLATEventCube::ebins(void) const
{
    return m_map.nmaps();
}


/***********************************************************************//**
 * @brief Return number of diffuse model components
 *
 * @return Number of diffuse model components.
 ***************************************************************************/
inline
int GLATEventCube::ndiffrsp(void) const
{
    return m_srcmap.size();
}

#endif /* GLATEVENTCUBE_HPP */
