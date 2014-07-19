/***************************************************************************
 *                  GCTAEventBin.hpp - CTA event bin class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2014 by Juergen Knoedlseder                         *
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
 * @file GCTAEventBin.hpp
 * @brief CTA event bin class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTAEVENTBIN_HPP
#define GCTAEVENTBIN_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GEventBin.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GCTAInstDir.hpp"


/***********************************************************************//**
 * @class GCTAEventBin
 *
 * @brief GCTAEventBin class interface definition
 *
 * This class implements a CTA event bin. The event bin is a collection of
 * pointers that points to the relevant information in memory. The class
 * itself does not allocate any memory, it just is a vector to collect
 * all relevant event information in a single place. This avoids duplication
 * of information.
 *
 * Setting up the pointers is done by the corresponding event bin container
 * class (GCTAEventCube).
 ***************************************************************************/
class GCTAEventBin : public GEventBin {

    // Friend classes
    friend class GCTAEventCube;

public:
    // Constructors and destructors
    GCTAEventBin(void);
    GCTAEventBin(const GCTAEventBin& bin);
    virtual ~GCTAEventBin(void);

    // Operators
    virtual GCTAEventBin& operator=(const GCTAEventBin& bin);

    // Implemented pure virtual base class methods
    virtual void               clear(void);
    virtual GCTAEventBin*      clone(void) const;
    virtual double             size(void) const;
    virtual const GCTAInstDir& dir(void) const;
    virtual const GEnergy&     energy(void) const;
    virtual const GTime&       time(void) const;
    virtual double             counts(void) const;
    virtual double             error(void) const;
    virtual void               counts(const double& counts);
    virtual std::string        print(const GChatter& chatter = NORMAL) const;

    // Other methods
    const double&  solidangle(void) const;
    const GEnergy& ewidth(void) const;
    const double&  ontime(void) const;
    const int&     ipix(void) const;
    const int&     ieng(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCTAEventBin& bin);
    void free_members(void);

    // Protected members
    int          m_ipix;        //!< Index in spatial map
    int          m_ieng;        //!< Index of energy layer
    GEnergy*     m_energy;      //!< Pointer to bin energy
    GCTAInstDir* m_dir;         //!< Pointer to bin direction
    GTime*       m_time;        //!< Pointer to bin time
    double*      m_counts;      //!< Pointer to number of counts
    double*      m_solidangle;  //!< Pointer to solid angle of pixel (sr)
    GEnergy*     m_ewidth;      //!< Pointer to energy width of bin
    double*      m_ontime;      //!< Pointer to ontime of bin (seconds)
};


/***********************************************************************//**
 * @brief Return the spatial pixel index
 *
 * @return Spatial pixel index.
 *
 * Returns the spatial index of the event bin in the event cube. If the event
 * is not part of an event cube, -1 is returned.
 ***************************************************************************/
inline
const int& GCTAEventBin::ipix(void) const
{
    return (m_ipix);
}


/***********************************************************************//**
 * @brief Return the energy layer index
 *
 * @return Energy layer index.
 *
 * Returns the energy layer of the event bin in the event cube. If the event
 * is not part of an event cube, -1 is returned.
 ***************************************************************************/
inline
const int& GCTAEventBin::ieng(void) const
{
    return (m_ieng);
}

#endif /* GCTAEVENTBIN_HPP */
