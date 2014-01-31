/***************************************************************************
 *             GCTAInstDir.hpp - CTA instrument direction class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2013 by Juergen Knoedlseder                         *
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
 * @file GCTAInstDir.hpp
 * @brief CTA instrument direction class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTAINSTDIR_HPP
#define GCTAINSTDIR_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GInstDir.hpp"
#include "GSkyDir.hpp"


/***********************************************************************//**
 * @class GCTAInstDir
 *
 * @brief CTA instrument direction class.
 *
 * The CTA instrument direction is an encapsulation of a sky direction
 * as CTA is an imaging device.
 ***************************************************************************/
class GCTAInstDir : public GInstDir {

public:
    // Constructors and destructors
    GCTAInstDir(void);
    explicit GCTAInstDir(const GSkyDir& dir);
    GCTAInstDir(const GCTAInstDir& dir);
    virtual ~GCTAInstDir(void);

    // Operators
    GCTAInstDir& operator=(const GCTAInstDir& dir);

    // Implemented pure virtual base class methods
    virtual void         clear(void);
    virtual GCTAInstDir* clone(void) const;
    virtual std::string  print(const GChatter& chatter = NORMAL) const;

    // Other methods
    void           dir(const GSkyDir& dir);
    GSkyDir&       dir(void);
    const GSkyDir& dir(void) const;
    void detx(const double &x){m_detx = x;};
    void dety(const double &y){m_dety = y;};
    const double& detx(void) const  {return m_detx;};
    const double& dety(void) const {return m_dety;};


protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCTAInstDir& dir);
    void free_members(void);

    // Data members
    GSkyDir m_dir;  //!< Observed incident direction of event
    double m_detx; //!< Instrument coordinate X
    double m_dety; //!< Instrument coordinate Y



};


/***********************************************************************//**
 * @brief Returns reference to sky direction
 *
 * @return Reference to sky direction.
 *
 * Returns reference to the sky direction.
 ***************************************************************************/
inline
GSkyDir& GCTAInstDir::dir(void)
{
    return (m_dir);
}


/***********************************************************************//**
 * @brief Returns reference to sky direction (const version)
 *
 * @return Reference to sky direction.
 *
 * Returns reference to the sky direction.
 ***************************************************************************/
inline
const GSkyDir& GCTAInstDir::dir(void) const
{
    return (m_dir);
}


/***********************************************************************//**
 * @brief Set sky direction
 *
 * @param[in] dir Sky direction.
 *
 * Set the sky direction.
 ***************************************************************************/
inline
void GCTAInstDir::dir(const GSkyDir& dir)
{
    m_dir = dir;
    return;
}

#endif /* GCTAINSTDIR_HPP */
