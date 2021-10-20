/***************************************************************************
 *               GXXXEventBin.hpp - [INSTRUMENT] event bin class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) [YEAR] by [AUTHOR]                                       *
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
 * @brief [INSTRUMENT] event bin class definition
 * @author [AUTHOR]
 */

#ifndef GXXXEVENTBIN_HPP
#define GXXXEVENTBIN_HPP

/* __ Includes ___________________________________________________________ */
#include "GEventBin.hpp"
#include "GXXXInstDir.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GPolarization.hpp"

/* __ Forward declarations _______________________________________________ */

/* __ Constants __________________________________________________________ */


/***********************************************************************//**
 * @class GXXXEventBin
 *
 * @brief [INSTRUMENT] event bin class
 *
 * This class defines an event bin of the [INSTRUMENT] event cube.
 *
 * Since many event bins share the same attributes (for example many bins
 * will actually have the same energy), it would be a waste of memory to
 * store all bin attributes together with the bin. Therefore, the
 * GXXXEventBin class implement it's own memory management. Either the 
 * class allocates memory for all attributes, or it takes pointers to
 * GXXXEventCube class member that store the information in an efficient
 * way.
 *
 * The data member m_alloc signals in which mode the bin operates. I true,
 * GXXXEventBin has allocated itself the memory for all attributes, and
 * hence has to take care about the memory allocation upon destruction.
 * Otherwise the pointers are just released.
 ***************************************************************************/
class GXXXEventBin : public GEventBin {

    // Friend classes
    friend class GXXXEventCube;

public:
    // Constructors and destructors
    GXXXEventBin(void);
    GXXXEventBin(const GXXXEventBin& bin);
    virtual ~GXXXEventBin(void);

    // Operators
    virtual GXXXEventBin& operator=(const GXXXEventBin& bin);

    // Implemented pure virtual base class methods
    virtual void                 clear(void);
    virtual GXXXEventBin*        clone(void) const;
    virtual std::string          classname(void) const;
    virtual double               size(void) const;
    virtual const GXXXInstDir&   dir(void) const;
    virtual const GEnergy&       energy(void) const;
    virtual const GTime&         time(void) const;
    virtual const GPolarization& polarization(void) const;
    virtual double               counts(void) const;
    virtual double               error(void) const;
    virtual void                 counts(const double& counts);
    virtual std::string          print(const GChatter& chatter = NORMAL) const;

    // Other methods
    // TODO: Add any further methods that are needed

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GXXXEventBin& bin);
    void free_members(void);

    // Protected members
    bool           m_alloc;        //!< Signals proper memory allocation
    int            m_index;        //!< Dataspace index
    double*        m_counts;       //!< Pointer to number of counts
    GXXXInstDir*   m_dir;          //!< Pointer to bin direction
    GEnergy*       m_energy;       //!< Pointer to bin energy
    GTime*         m_time;         //!< Pointer to bin time
    GPolarization* m_polarization; //!< Pointer to bin polarization
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GXXXEventBin").
 ***************************************************************************/
inline
std::string GXXXEventBin::classname(void) const
{
    return ("GXXXEventBin");
}


/***********************************************************************//**
 * @brief Return instrument direction
 *
 * @return Instrument direction.
 *
 * Returns the instrument direction of the event bin.
 ***************************************************************************/
inline
const GXXXInstDir& GXXXEventBin::dir(void) const
{
    return (*m_dir);
}


/***********************************************************************//**
 * @brief Return energy
 *
 * @return Energy.
 *
 * Returns the energy of the event bin.
 ***************************************************************************/
inline
const GEnergy& GXXXEventBin::energy(void) const
{
    return (*m_energy);
}


/***********************************************************************//**
 * @brief Return time
 *
 * @return Time.
 *
 * Returns the time of the event bin.
 ***************************************************************************/
inline
const GTime& GXXXEventBin::time(void) const
{
    return (*m_time);
}


/***********************************************************************//**
 * @brief Return polarization
 *
 * @return Polarization.
 *
 * Returns the polarization of the event bin.
 ***************************************************************************/
inline
const GPolarization& GXXXEventBin::polarization(void) const
{
    return (*m_polarization);
}


/***********************************************************************//**
 * @brief Return number of counts
 *
 * @return Number of counts.
 *
 * Returns the number of counts in the event bin.
 ***************************************************************************/
inline
double GXXXEventBin::counts(void) const
{
    return (*m_counts);
}


/***********************************************************************//**
 * @brief Set number of counts
 *
 * @param[in] counts Number of counts.
 *
 * Set the number of counts in the event bin.
 ***************************************************************************/
inline
void GXXXEventBin::counts(const double& counts)
{
    *m_counts = counts;
    return;
}

#endif /* GXXXEVENTBIN_HPP */
