/***************************************************************************
 *             GCTAEventCube.hpp - CTA event bin container class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2016 by Juergen Knoedlseder                         *
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
#include "GSkyMap.hpp"
#include "GEbounds.hpp"
#include "GGti.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GCTAInstDir.hpp"
#include "GFitsTable.hpp"
#include "GFitsImage.hpp"

/* __ Forward declarations _______________________________________________ */
class GFilename;


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
    explicit GCTAEventCube(const GFilename& filename);
    GCTAEventCube(const GSkyMap& map, const GEbounds& ebds, const GGti& gti);
    GCTAEventCube(const GSkyMap& map, const GSkyMap& weights,
                  const GEbounds& ebds, const GGti& gti);
    GCTAEventCube(const GCTAEventCube& cube);
    virtual ~GCTAEventCube(void);

    // Operators
    virtual GCTAEventCube&      operator=(const GCTAEventCube& cube);
    virtual GCTAEventBin*       operator[](const int& index);
    virtual const GCTAEventBin* operator[](const int& index) const;

    // Implemented pure virtual base class methods
    virtual void           clear(void);
    virtual GCTAEventCube* clone(void) const;
    virtual std::string    classname(void) const;
    virtual int            size(void) const;
    virtual int            dim(void) const;
    virtual int            naxis(const int& axis) const;
    virtual void           load(const GFilename& filename);
    virtual void           save(const GFilename& filename,
                                const bool&      clobber = false) const;
    virtual void           read(const GFits& file);
    virtual void           write(GFits& file) const;
    virtual int            number(void) const;
    virtual std::string    print(const GChatter& chatter = NORMAL) const;

    // Other methods
    const GTime&           time(void) const;
    const GEnergy&         energy(const int& index) const;
    void                   counts(const GSkyMap& counts);
    const GSkyMap&         counts(void) const;
    void                   weights(const GSkyMap& weights);
    const GSkyMap&         weights(void) const;
    int                    nx(void) const;
    int                    ny(void) const;
    int                    npix(void) const;
    int                    ebins(void) const;

protected:
    // Protected methods
    void         init_members(void);
    void         copy_members(const GCTAEventCube& cube);
    void         free_members(void);
    void         read_cntmap(const GFitsImage& hdu);
    void         read_ebds(const GFitsTable& hdu);
    void         read_gti(const GFitsTable& hdu);
    void         set_directions(void);
    virtual void set_energies(void);
    virtual void set_times(void);
    void         init_bin(void);
    void         set_bin(const int& index);

    // Protected members
    GSkyMap                  m_map;        //!< Counts cube stored as sky map
    GSkyMap                  m_weights;    //!< Cube weights stored as sky map
    GCTAEventBin             m_bin;        //!< Actual event bin
    GTime                    m_time;       //!< Event cube mean time
    std::vector<GCTAInstDir> m_dirs;       //!< Array of event directions
    std::vector<double>      m_solidangle; //!< Array of solid angles (sr)
    std::vector<GEnergy>     m_energies;   //!< Array of log mean energies
    std::vector<GEnergy>     m_ewidth;     //!< Array of energy bin widths
    double                   m_ontime;     //!< Event cube ontime (sec)
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCTAEventCube").
 ***************************************************************************/
inline
std::string GCTAEventCube::classname(void) const
{
    return ("GCTAEventCube");
}


/***********************************************************************//**
 * @brief Return dimension of event cube
 *
 * @return Dimension of event cube.
 ***************************************************************************/
inline
int GCTAEventCube::dim(void) const
{
    // Return dimension
    return ((m_map.nmaps() > 0) ? 3 : 0);
}


/***********************************************************************//**
 * @brief Return event cube counts as sky map
 *
 * @return Event cube counts.
 *
 * Returns the event cube counts as sky map.
 ***************************************************************************/
inline
const GSkyMap& GCTAEventCube::counts(void) const
{
    return (m_map);
}


/***********************************************************************//**
 * @brief Return event cube weights as sky map
 *
 * @return Event cube weights.
 *
 * Returns the event cube weights as sky map.
 ***************************************************************************/
inline
const GSkyMap& GCTAEventCube::weights(void) const
{
    return (m_weights);
}


/***********************************************************************//**
 * @brief Set event cube weights from sky map
 *
 * @param[in] weights Event cube weights sky map.
 *
 * Sets the event cube weights from sky map.
 ***************************************************************************/
inline
void GCTAEventCube::weights(const GSkyMap& weights)
{
    m_weights = weights;
    return;
}


/***********************************************************************//**
 * @brief Return number of bins in X direction
 *
 * @return Number of bins in X direction.
 ***************************************************************************/
inline
int GCTAEventCube::nx(void) const
{
    return (m_map.nx());
}


/***********************************************************************//**
 * @brief Return number of bins in Y direction
 *
 * @return Number of bins in Y direction.
 ***************************************************************************/
inline
int GCTAEventCube::ny(void) const
{
    return (m_map.ny());
}


/***********************************************************************//**
 * @brief Return number of pixels in one energy bins of the event cube
 *
 * @return Number of pixels in one energy bins of the event cube.
 ***************************************************************************/
inline
int GCTAEventCube::npix(void) const
{
    return (m_map.npix());
}


/***********************************************************************//**
 * @brief Return number of energy bins in the event cube
 *
 * @return Number of energy bins in the event cube.
 ***************************************************************************/
inline
int GCTAEventCube::ebins(void) const
{
    return (m_map.nmaps());
}


/***********************************************************************//**
 * @brief Return event cube mean time
 *
 * @return Event cube mean time.
 ***************************************************************************/
inline
const GTime& GCTAEventCube::time(void) const
{
    return (m_time);
}

#endif /* GCTAEVENTCUBE_HPP */
