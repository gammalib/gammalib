/***************************************************************************
 *                 GHorizDir.hpp - Horizontal direction class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Karl Kosack                                      *
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
 * @file GHorizDir.hpp
 * @brief Horizontal direction class interface definition
 * @author Karl Kosack 
 */

#ifndef GHORIZDIR_HPP
#define GHORIZDIR_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GVector.hpp"
#include "GMath.hpp"

/* __ Compile options ____________________________________________________ */


/***********************************************************************//**
 * @class GHorizDir
 *
 * @brief Horizontal (Alt/Az) direction class
 *
 * This class is essentially a copy of GSkyDir and implements a
 * spherical coordinate on the sky, in horizontal coordinates as seen
 * from Earth. 
 *
 * @note since this class shares much functionality with GSkyDir, a
 * common base-class should probably be created in the future to avoid
 * duplication of code.
 ***************************************************************************/
class GHorizDir : public GBase {

    // Operator friends
    friend bool operator==(const GHorizDir &a, const GHorizDir &b);
    friend bool operator!=(const GHorizDir &a, const GHorizDir &b);

public:
    // Constructors and destructors
    GHorizDir(void);
    GHorizDir(const GHorizDir& dir);
    virtual ~GHorizDir(void);

    // Operators
    GHorizDir& operator=(const GHorizDir& dir);

    // Methods
    void          clear(void);
    GHorizDir*    clone(void) const;
    std::string   classname(void) const;
    void          altaz(const double& alt, const double& az);
    void          altaz_deg(const double& alt, const double& az);
    void          celvector(const GVector& vector);
    void          rotate_deg(const double& phi, const double& theta);
    const double& alt(void) const;
    const double& az(void) const;
    double        zenith(void) const;
    double        zenith_deg(void) const;
    double        alt_deg(void) const;
    double        az_deg(void) const;
    GVector       celvector(void) const;
    double        dist(const GHorizDir& dir) const;
    double        dist_deg(const GHorizDir& dir) const;
    std::string   print(const GChatter& chatter = NORMAL) const;

private:
    // Private methods
    void init_members(void);
    void copy_members(const GHorizDir& dir);
    void free_members(void);

    // Private members
    double m_alt;       //!< altitude in radians
    double m_az;        //!< azimuth in radians
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GHorizDir").
 ***************************************************************************/
inline
std::string GHorizDir::classname(void) const
{
    return ("GHorizDir");
}


/***********************************************************************//**
 * @brief Return zenith angle in radians
 *
 * @return zenith angle in radians.
 ***************************************************************************/
inline
double GHorizDir::zenith() const
{
    return (gammalib::pihalf - m_alt);
}


/***********************************************************************//**
 * @brief Return zenith angle in degrees
 *
 * @return zenith angle in degrees.
 ***************************************************************************/
inline
double GHorizDir::zenith_deg() const
{
    return (zenith() * gammalib::rad2deg);
}


/***********************************************************************//**
 * @brief Return altitude angle in radians
 *
 * @return Altitude angle in radians.
 ***************************************************************************/
inline
const double& GHorizDir::alt() const
{
    return m_alt;
}


/***********************************************************************//**
 * @brief Return altitude angle in degrees
 *
 * @return Altitude angle in degrees.
 ***************************************************************************/
inline
double GHorizDir::alt_deg() const 
{
    return (m_alt * gammalib::rad2deg);
}


/***********************************************************************//**
 * @brief Return azimuth angle in radians
 *
 * @return Azimuth angle in radians.
 ***************************************************************************/
inline
const double& GHorizDir::az() const
{
    return m_az;
}


/***********************************************************************//**
 * @brief Return azimuth angle in degrees
 *
 * @return Azimuth angle in degrees.
 ***************************************************************************/
inline
double GHorizDir::az_deg() const 
{
    return (m_az * gammalib::rad2deg);
}

#endif /* GHORIZDIR_HPP */
