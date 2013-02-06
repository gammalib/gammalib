/***************************************************************************
 *                 GSPIEventBin.hpp  -  SPI event bin class                *
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
 * @file GSPIEventBin.hpp
 * @brief SPI event bin class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GSPIEVENTBIN_HPP
#define GSPIEVENTBIN_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GEventBin.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GSPIInstDir.hpp"


/***********************************************************************//**
 * @class GSPIEventBin
 *
 * @brief SPI event bin class
 ***************************************************************************/
class GSPIEventBin : public GEventBin {

    // Friend classes
    //friend class GSPIEventCube;

public:
    // Constructors and destructors
    GSPIEventBin(void);
    GSPIEventBin(const GSPIEventBin& bin);
    virtual ~GSPIEventBin(void);

    // Operators
    virtual GSPIEventBin& operator=(const GSPIEventBin& bin);

    // Implemented pure virtual base class methods
    virtual void               clear(void);
    virtual GSPIEventBin*      clone(void) const;
    virtual double             size(void) const;
    virtual const GSPIInstDir& dir(void) const { return *m_dir; }
    virtual const GEnergy&     energy(void) const { return *m_energy; }
    virtual const GTime&       time(void) const { return *m_time; }
    virtual double             counts(void) const { return *m_counts; }
    virtual double             error(void) const;
    virtual void               counts(const double& counts) { *m_counts=counts; }
    virtual std::string        print(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GSPIEventBin& bin);
    void free_members(void);

    // Protected members
    bool         m_alloc;       //!< Signals proper memory allocation
    int          m_index;       //!< Dataspace index
    double*      m_counts;      //!< Pointer to number of counts
    GSPIInstDir* m_dir;         //!< Pointer to bin direction
    GEnergy*     m_energy;      //!< Pointer to bin energy
    GTime*       m_time;        //!< Pointer to bin time
};

#endif /* GSPIEVENTBIN_HPP */
