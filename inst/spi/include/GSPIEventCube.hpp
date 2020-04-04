/***************************************************************************
 *            GSPIEventCube.hpp - INTEGRAL/SPI event cube class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2020 by Juergen Knoedlseder                              *
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
 * @file GSPIEventCube.hpp
 * @brief INTEGRAL/SPI event bin container class definition
 * @author Juergen Knoedlseder
 */

#ifndef GSPIEVENTCUBE_HPP
#define GSPIEVENTCUBE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GEventCube.hpp"
#include "GSPIEventBin.hpp"
#include "GSPIInstDir.hpp"
#include "GSkyDir.hpp"
#include "GSkyMap.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"

/* __ Forward declarations _______________________________________________ */
class GFilename;
class GEbounds;
class GGti;
class GFits;

/* __ Constants __________________________________________________________ */


/***********************************************************************//**
 * @class GSPIEventCube
 *
 * @brief INTEGRAL/SPI event bin container class
 *
 * This class is a container class for INTEGRAL/SPI event bins.
 ***************************************************************************/
class GSPIEventCube : public GEventCube {

public:
    // Constructors and destructors
    GSPIEventCube(void);
    explicit GSPIEventCube(const GFilename& filename);
    GSPIEventCube(const GSPIEventCube& cube);
    virtual ~GSPIEventCube(void);

    // Operators
    virtual GSPIEventCube&      operator=(const GSPIEventCube& cube);
    virtual GSPIEventBin*       operator[](const int& index);
    virtual const GSPIEventBin* operator[](const int& index) const;

    // Implemented pure virtual base class methods
    virtual void           clear(void);
    virtual GSPIEventCube* clone(void) const;
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
    // TODO: Add any further methods that are needed

protected:
    // Protected methods
    void         init_members(void);
    void         copy_members(const GSPIEventCube& cube);
    void         free_members(void);
    virtual void set_energies(void);
    virtual void set_times(void);
    void         init_bin(void);
    void         set_bin(const int& index);

    // Protected members
    GSPIEventBin m_bin;        //!< Actual event bin
    GSPIInstDir  m_dir;        //!< Actual event direction
    GTime        m_time;       //!< Event cube mean time
    double       m_ontime;     //!< Event cube ontime (sec)
    GEnergy      m_energy;     //!< Event cube mean energy
    GEnergy      m_ewidth;     //!< Event cube energy bin width
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GSPIEventCube").
 ***************************************************************************/
inline
std::string GSPIEventCube::classname(void) const
{
    return ("GSPIEventCube");
}


/***********************************************************************//**
 * @brief Set energies
 *
 * Sets energies.
 ***************************************************************************/
inline
void GSPIEventCube::set_energies(void)
{
    return;
}


/***********************************************************************//**
 * @brief Set times
 *
 * Sets times.
 ***************************************************************************/
inline
void GSPIEventCube::set_times(void)
{
    return;
}

#endif /* GSPIEVENTCUBE_HPP */
