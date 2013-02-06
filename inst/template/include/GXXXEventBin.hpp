/***************************************************************************
 *                 GXXXEventBin.hpp  -  XXX event bin class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013 by Juergen Knoedlseder                              *
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
 * @file GXXXEventBin.hpp
 * @brief XXX event bin class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GXXXEVENTBIN_HPP
#define GXXXEVENTBIN_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GEventBin.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GXXXInstDir.hpp"


/***********************************************************************//**
 * @class GXXXEventBin
 *
 * @brief XXX event bin class
 *
 * This class defines an event bin of the XXX data space. The class holds
 * pointers to the bin attributes, and allocates memory so that all bin
 * attributes have an associated memory space. With this technique, the class
 * can be used to hold information that is external to GXXXEventBin. The
 * GXXXEventCube class will use this property to directly manipulate the
 * pointers of the bin.
 ***************************************************************************/
class GXXXEventBin : public GEventBin {

public:
    // Constructors and destructors
    GXXXEventBin(void);
    GXXXEventBin(const GXXXEventBin& bin);
    virtual ~GXXXEventBin(void);

    // Operators
    virtual GXXXEventBin& operator=(const GXXXEventBin& bin);

    // Implemented pure virtual base class methods
    virtual void               clear(void);
    virtual GXXXEventBin*      clone(void) const;
    virtual double             size(void) const;
    virtual const GXXXInstDir& dir(void) const { return *m_dir; }
    virtual const GEnergy&     energy(void) const { return *m_energy; }
    virtual const GTime&       time(void) const { return *m_time; }
    virtual double             counts(void) const { return *m_counts; }
    virtual double             error(void) const;
    virtual void               counts(const double& counts) { *m_counts=counts; }
    virtual std::string        print(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GXXXEventBin& bin);
    void free_members(void);

    // Protected members
    bool         m_alloc;    //!< Signals memory allocation by object
    int          m_index;    //!< Dataspace index (-1 if not used)
    double*      m_counts;   //!< Pointer to number of counts
    GXXXInstDir* m_dir;      //!< Pointer to bin direction
    GEnergy*     m_energy;   //!< Pointer to bin energy
    GTime*       m_time;     //!< Pointer to bin time
};

#endif /* GXXXEVENTBIN_HPP */
