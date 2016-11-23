/***************************************************************************
 *                GCOMEventBin.hpp - COMPTEL event bin class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2016 by Juergen Knoedlseder                         *
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
 * @file GCOMEventBin.hpp
 * @brief COMPTEL event bin class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GCOMEVENTBIN_HPP
#define GCOMEVENTBIN_HPP

/* __ Includes ___________________________________________________________ */
#include "GEventBin.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GCOMInstDir.hpp"


/***********************************************************************//**
 * @class GCOMEventBin
 *
 * @brief COMPTEL event bin class
 *
 * This class defines an event bin of the COMPTEL data space. The class holds
 * pointers to the bin attributes, and allocates memory so that all bin
 * attributes have an associated memory space. With this technique, the class
 * can be used to hold information that is external to GCOMEventBin. The
 * GCOMEventCube class will use this property to directly manipulate the
 * pointers of the bin.
 ***************************************************************************/
class GCOMEventBin : public GEventBin {

    // Friend classes
    friend class GCOMEventCube;

public:
    // Constructors and destructors
    GCOMEventBin(void);
    GCOMEventBin(const GCOMEventBin& bin);
    virtual ~GCOMEventBin(void);

    // Operators
    virtual GCOMEventBin& operator= (const GCOMEventBin& bin);

    // Implemented pure virtual base class methods
    virtual void               clear(void);
    virtual GCOMEventBin*      clone(void) const;
    virtual std::string        classname(void) const;
    virtual double             size(void) const;
    virtual const GCOMInstDir& dir(void) const;
    virtual const GEnergy&     energy(void) const;
    virtual const GTime&       time(void) const;
    virtual double             counts(void) const;
    virtual double             error(void) const;
    virtual void               counts(const double& counts);
    virtual std::string        print(const GChatter& chatter = NORMAL) const;

    // Other methods
    const int&     index(void) const;
    const double&  solidangle(void) const;
    const GEnergy& ewidth(void) const;
    const double&  ontime(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCOMEventBin& bin);
    void free_members(void);

    // Protected members
    bool         m_alloc;       //!< Signals proper memory allocation
    int          m_index;       //!< Dataspace index
    GCOMInstDir* m_dir;         //!< Pointer to bin direction
    GTime*       m_time;        //!< Pointer to bin time
    GEnergy*     m_energy;      //!< Pointer to bin energy
    GEnergy*     m_ewidth;      //!< Pointer to energy width of bin
    double*      m_counts;      //!< Pointer to number of counts
    double*      m_solidangle;  //!< Pointer to solid angle of pixel (sr)
    double*      m_ontime;      //!< Pointer to ontime of bin (seconds)
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCOMEventBin").
 ***************************************************************************/
inline
std::string GCOMEventBin::classname(void) const
{
    return ("GCOMEventBin");
}


/***********************************************************************//**
 * @brief Return instrument direction
 *
 * @return Instrument direction.
 *
 * Returns the instrument direction of the event bin.
 ***************************************************************************/
inline
const GCOMInstDir& GCOMEventBin::dir(void) const
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
const GEnergy& GCOMEventBin::energy(void) const
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
const GTime& GCOMEventBin::time(void) const
{
    return (*m_time);
}


/***********************************************************************//**
 * @brief Return number of counts
 *
 * @return Number of counts.
 *
 * Returns the number of counts in the event bin.
 ***************************************************************************/
inline
double GCOMEventBin::counts(void) const
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
void GCOMEventBin::counts(const double& counts)
{
    *m_counts = counts;
    return;
}


/***********************************************************************//**
 * @brief Return bin index
 *
 * @return Bin index.
 *
 * Returns the index of the event bin if the bin is part of an event cube. If
 * the event is not part of an event cube, -1 is returned.
 ***************************************************************************/
inline
const int& GCOMEventBin::index(void) const
{
    return (m_index);
}


/***********************************************************************//**
 * @brief Return solid angle
 *
 * @return Solid angle.
 *
 * Returns the solid angle of the event bin.
 ***************************************************************************/
inline
const double& GCOMEventBin::solidangle(void) const
{
    return (*m_solidangle);
}


/***********************************************************************//**
 * @brief Return energy width
 *
 * @return Energy width.
 *
 * Returns the energy width of the event bin.
 ***************************************************************************/
inline
const GEnergy& GCOMEventBin::ewidth(void) const
{
    return (*m_ewidth);
}


/***********************************************************************//**
 * @brief Return ontime
 *
 * @return Ontime.
 *
 * Returns the ontime of the event bin.
 ***************************************************************************/
inline
const double& GCOMEventBin::ontime(void) const
{
    return (*m_ontime);
}

#endif /* GCOMEVENTBIN_HPP */
