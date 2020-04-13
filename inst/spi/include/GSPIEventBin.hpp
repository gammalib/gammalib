/***************************************************************************
 *             GSPIEventBin.hpp - INTEGRAL/SPI event bin class             *
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
 * @file GSPIEventBin.hpp
 * @brief INTEGRAL/SPI event bin class definition
 * @author Juergen Knoedlseder
 */

#ifndef GSPIEVENTBIN_HPP
#define GSPIEVENTBIN_HPP

/* __ Includes ___________________________________________________________ */
#include "GEventBin.hpp"
#include "GSPIInstDir.hpp"

/* __ Forward declarations _______________________________________________ */
class GTime;
class GEnergy;
class GSPIEventCube;

/* __ Constants __________________________________________________________ */


/***********************************************************************//**
 * @class GSPIEventBin
 *
 * @brief INTEGRAL/SPI event bin class
 *
 * This class defines an event bin of the INTEGRAL/SPI event cube.
 *
 * Since many event bins share the same attributes (for example many bins
 * will actually have the same energy), it would be a waste of memory to
 * store all bin attributes together with the bin. Therefore, the
 * GSPIEventBin class implement it's own memory management. Either the 
 * class allocates memory for all attributes, or it takes pointers to
 * GSPIEventCube class member that store the information in an efficient
 * way.
 *
 * The data member m_alloc signals in which mode the bin operates. I true,
 * GSPIEventBin has allocated itself the memory for all attributes, and
 * hence has to take care about the memory allocation upon destruction.
 * Otherwise the pointers are just released.
 ***************************************************************************/
class GSPIEventBin : public GEventBin {

    // Friend classes
    friend class GSPIEventCube;

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
    virtual std::string        classname(void) const;
    virtual double             size(void) const;
    virtual const GSPIInstDir& dir(void) const;
    virtual const GEnergy&     energy(void) const;
    virtual const GTime&       time(void) const;
    virtual double             counts(void) const;
    virtual double             error(void) const;
    virtual void               counts(const double& counts);
    virtual std::string        print(const GChatter& chatter = NORMAL) const;

    // Other methods
    const double& model(const int& index) const;
    const double& ontime(void) const;
    const double& livetime(void) const;
    const int&    index(void) const;
    const int&    ipt(void) const;
    const int&    idir(void) const;
    const int&    iebin(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GSPIEventBin& bin);
    void free_members(void);

    // Protected members
    bool         m_alloc;      //!< Signals proper memory allocation
    int          m_index;      //!< Dataspace index
    int          m_ipt;        //!< Pointing index
    int          m_idir;       //!< Direction index
    int          m_iebin;      //!< Energy bin index
    int          m_num_models; //!< Number of models in bin
    GSPIInstDir* m_dir;        //!< Pointer to direction of bin
    GTime*       m_time;       //!< Pointer to time of bin
    GEnergy*     m_energy;     //!< Pointer to energy of bin
    double*      m_counts;     //!< Pointer to number of counts
    double*      m_ontime;     //!< Pointer to ontime of bin
    double*      m_livetime;   //!< Pointer to livetime of bin
    double*      m_size;       //!< Pointer to size of bin
    double*      m_models;     //!< Pointer to models of bin
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GSPIEventBin").
 ***************************************************************************/
inline
std::string GSPIEventBin::classname(void) const
{
    return ("GSPIEventBin");
}


/***********************************************************************//**
 * @brief Return instrument direction
 *
 * @return Instrument direction.
 *
 * Returns the instrument direction of the event bin.
 ***************************************************************************/
inline
const GSPIInstDir& GSPIEventBin::dir(void) const
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
const GEnergy& GSPIEventBin::energy(void) const
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
const GTime& GSPIEventBin::time(void) const
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
double GSPIEventBin::counts(void) const
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
void GSPIEventBin::counts(const double& counts)
{
    *m_counts = counts;
    return;
}

/***********************************************************************//**
 * @brief Return size of event bin
 *
 * @return Size of event bin (MeV s)
 *
 * Returns the size of the event bin.
 ***************************************************************************/
inline
double GSPIEventBin::size(void) const
{
    return (*m_size);
}


/***********************************************************************//**
 * @brief Return ontime of event bin
 *
 * @return Size of ontime of bin (s)
 *
 * Returns the ontime of the event bin.
 ***************************************************************************/
inline
const double& GSPIEventBin::ontime(void) const
{
    return (*m_ontime);
}


/***********************************************************************//**
 * @brief Return livetime of event bin
 *
 * @return Size of livetime of bin (s)
 *
 * Returns the livetime of the event bin.
 ***************************************************************************/
inline
const double& GSPIEventBin::livetime(void) const
{
    return (*m_livetime);
}


/***********************************************************************//**
 * @brief Return event bin index
 *
 * @return Event bin index
 *
 * Returns the event bin index in the event cube.
 ***************************************************************************/
inline
const int& GSPIEventBin::index(void) const
{
    return (m_index);
}


/***********************************************************************//**
 * @brief Return event bin pointing index
 *
 * @return Event bin pointing index
 *
 * Returns the event bin pointing index in the event cube.
 ***************************************************************************/
inline
const int& GSPIEventBin::ipt(void) const
{
    return (m_ipt);
}


/***********************************************************************//**
 * @brief Return event bin direction index
 *
 * @return Event bin direction index
 *
 * Returns the event bin instrument direction index in the event cube.
 ***************************************************************************/
inline
const int& GSPIEventBin::idir(void) const
{
    return (m_idir);
}


/***********************************************************************//**
 * @brief Return event bin energy index
 *
 * @return Event bin energy index
 *
 * Returns the event bin energy index in the event cube.
 ***************************************************************************/
inline
const int& GSPIEventBin::iebin(void) const
{
    return (m_iebin);
}

#endif /* GSPIEVENTBIN_HPP */
