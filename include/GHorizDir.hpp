/***************************************************************************
 *                  GHorizDir.hpp - Horizontal direction class                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2013 by Juergen Knoedlseder                         *
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
 * @author K. Kosack
 */

#ifndef GSKYDIR_HPP
#define GSKYDIR_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GVector.hpp"

/* __ Compile options ____________________________________________________ */
#define G_SINCOS_CACHE


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
    GHorizDir& operator= (const GHorizDir& dir);

    // Methods
    void          clear(void);
    GHorizDir*    clone(void) const;
    void          altaz(const double& alt, const double& az);
    void          altaz_deg(const double& alt, const double& az);
    void          celvector(const GVector& vector);
    void          rotate_deg(const double& phi, const double& theta);
    const double& alt(void) const;
    const double& az(void) const;
    double        zenith() const;
    double        zenith_deg() const;
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
    double m_alt;         //!< altitude in radians
    double m_az;        //!< azimuth in radians

    // Sincos cache
    // #if defined(G_SINCOS_CACHE)
    // mutable bool   m_has_lb_cache;
    // mutable bool   m_has_radec_cache;
    // mutable double m_sin_b;
    // mutable double m_cos_b;
    // mutable double m_sin_dec;
    // mutable double m_cos_dec;
    // #endif
};


/***********************************************************************//**
 * @brief Return zenith angle in radians
 *
 * @return zenith angle in radians
 ***************************************************************************/
inline
double GHorizDir::zenith() 
{
    return gammalib::pihalf - m_alt;
};

/***********************************************************************//**
 * @brief Return zenith angle in radians
 *
 * @return zenith angle in radians
 ***************************************************************************/
inline
double GHorizDir::zenith_deg() 
{
    return zenith() * gammalib::rad2deg;
};


#endif /* GHORIZDIR_HPP */
