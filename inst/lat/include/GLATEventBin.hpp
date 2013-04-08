/***************************************************************************
 *                 GLATEventBin.hpp - LAT event bin class                  *
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
 * @file GLATEventBin.hpp
 * @brief Fermi-LAT event bin class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GLATEVENTBIN_HPP
#define GLATEVENTBIN_HPP

/* __ Includes ___________________________________________________________ */
#include "GEventBin.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GLATInstDir.hpp"

/* __ Forward declarations _______________________________________________ */
class GLATEventCube;


/***********************************************************************//**
 * @class GLATEventBin
 *
 * @brief Fermi-LAT event bin class
 *
 * This class implement a counts map bin for the Fermi-LAT telescope.
 ***************************************************************************/
class GLATEventBin : public GEventBin {

    // Friend classes
    friend class GLATEventCube;

public:
    // Constructors and destructors
    GLATEventBin(void);
    GLATEventBin(const GLATEventBin& bin);
    virtual ~GLATEventBin(void);

    // Operators
    virtual GLATEventBin& operator= (const GLATEventBin& bin);

    // Implemented pure virtual base class methods
    virtual void               clear(void);
    virtual GLATEventBin*      clone(void) const;
    virtual double             size(void) const;
    virtual const GLATInstDir& dir(void) const;
    virtual const GEnergy&     energy(void) const;
    virtual const GTime&       time(void) const;
    virtual double             counts(void) const;
    virtual double             error(void) const;
    virtual void               counts(const double& counts);
    virtual std::string        print(const GChatter& chatter = NORMAL) const;

    // Other methods
    const double&  omega(void) const;
    const GEnergy& ewidth(void) const;
    const double&  ontime(void) const;
    const int&     index(void) const { return m_index; }
    const int&     ipix(void) const { return m_ipix; }
    const int&     ieng(void) const { return m_ieng; }
    GLATEventCube* cube(void) const { return m_cube; }

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GLATEventBin& bin);
    void free_members(void);

    // Protected members
    GLATEventCube* m_cube;        //!< Event cube back pointer
    int            m_index;       //!< Actual skymap index
    int            m_ipix;        //!< Actual spatial index
    int            m_ieng;        //!< Actual energy index
    GEnergy*       m_energy;      //!< Pointer to bin energy
    GLATInstDir*   m_dir;         //!< Pointer to bin direction
    GTime*         m_time;        //!< Pointer to bin time
    double*        m_counts;      //!< Pointer to number of counts
    double*        m_omega;       //!< Pointer to solid angle of pixel (sr)
    GEnergy*       m_ewidth;      //!< Pointer to energy width of bin
    double*        m_ontime;      //!< Pointer to ontime of bin (seconds)
};

#endif /* GLATEVENTBIN_HPP */
