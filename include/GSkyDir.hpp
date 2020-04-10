/***************************************************************************
 *                     GSkyDir.hpp - Sky direction class                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2020 by Juergen Knoedlseder                         *
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
 * @file GSkyDir.hpp
 * @brief Sky direction class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GSKYDIR_HPP
#define GSKYDIR_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GVector.hpp"
#include "GMath.hpp"

/* __ Forward declarations _______________________________________________ */
class GTime;

/* __ Compile options ____________________________________________________ */
#define G_SINCOS_CACHE

/* __ Kluge for compilation on Sun _______________________________________ */
#ifdef sun
#undef sun
#endif


/***********************************************************************//**
 * @class GSkyDir
 *
 * @brief Sky direction class
 *
 * This class implements a spherical coordinate on the sky. Two coordinate
 * systems are supported: celestial (or equatorial), defined by Right
 * Ascension and Declination, and Galactic, defined by galactic longitude
 * and galactic latitude. The class provides automatic conversion between
 * both systems if required. Coordinates are stored in either of the
 * systems (in units of radians), and conversion is performed (and stored)
 * if requested. Coordinates can be given and returned in radians or in
 * degrees. Note that the epoch for celestial coordinates is fixed to J2000.
 ***************************************************************************/
class GSkyDir : public GBase {

    // Operator friends
    friend bool operator==(const GSkyDir &a, const GSkyDir &b);
    friend bool operator!=(const GSkyDir &a, const GSkyDir &b);

public:
    // Constructors and destructors
    GSkyDir(void);
    explicit GSkyDir(const GVector& vector);
    GSkyDir(const GSkyDir& dir);
    virtual ~GSkyDir(void);

    // Operators
    GSkyDir& operator=(const GSkyDir& dir);

    // Methods
    void          clear(void);
    GSkyDir*      clone(void) const;
    std::string   classname(void) const;
    void          radec(const double& ra, const double& dec);
    void          radec_deg(const double& ra, const double& dec);
    void          lb(const double& l, const double& b);
    void          lb_deg(const double& l, const double& b);
    void          celvector(const GVector& vector);
    void          rotate(const double& phi, const double& theta);
    void          rotate_deg(const double& phi, const double& theta);
    void          precess(const double& from_epoch, const double& to_epoch);
    void          sun(const GTime& time, const double& epoch = 2000.0);
    void          moon(const GTime& time, const double& epoch = 2000.0);
    const double& l(void) const;
    const double& b(void) const;
    const double& ra(void) const;
    const double& dec(void) const;
    double        l_deg(void) const;
    double        b_deg(void) const;
    double        ra_deg(void) const;
    double        dec_deg(void) const;
    GVector       celvector(void) const;
    double        cos_dist(const GSkyDir& dir) const;
    double        dist(const GSkyDir& dir) const;
    double        dist_deg(const GSkyDir& dir) const;
    double        posang(const GSkyDir& dir) const;
    double        posang_deg(const GSkyDir& dir) const;
    std::string   print(const GChatter& chatter = NORMAL) const;

private:
    // Private methods
    void init_members(void);
    void copy_members(const GSkyDir& dir);
    void free_members(void);
    void equ2gal(void) const;
    void gal2equ(void) const;
    void euler(const int& type, const double& xin, const double &yin,
               double* xout, double *yout) const;

    // Private members
    mutable bool   m_has_lb;     //!< Has galactic coordinates
    mutable bool   m_has_radec;  //!< Has equatorial coordinates
    double         m_l;          //!< Galactic longitude in radians
    double         m_b;          //!< Galactic latitude in radians
    double         m_ra;         //!< Right Ascension in radians
    double         m_dec;        //!< Declination in radians

    // Sincos cache
    #if defined(G_SINCOS_CACHE)
    mutable bool   m_has_lb_cache;
    mutable bool   m_has_radec_cache;
    mutable double m_sin_b;
    mutable double m_cos_b;
    mutable double m_sin_dec;
    mutable double m_cos_dec;
    #endif
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GSkyDir").
 ***************************************************************************/
inline
std::string GSkyDir::classname(void) const
{
    return ("GSkyDir");
}


/***********************************************************************//**
 * @brief Return galactic longitude in radians
 *
 * @return Galactic longitude in radians.
 ***************************************************************************/
inline
const double& GSkyDir::l(void) const
{
    if (!m_has_lb && m_has_radec) {
        equ2gal();
    }
    return m_l;
}


/***********************************************************************//**
 * @brief Return galactic longitude in degrees
 *
 * @return Galactic longitude in degrees.
 ***************************************************************************/
inline
double GSkyDir::l_deg(void) const
{
    return (l() * gammalib::rad2deg);
}


/***********************************************************************//**
 * @brief Return galactic latitude in radians
 *
 * @return Galactic latitude in radians.
 ***************************************************************************/
inline
const double& GSkyDir::b(void) const
{
    if (!m_has_lb && m_has_radec) {
        equ2gal();
    }
    return m_b;
}


/***********************************************************************//**
 * @brief Returns galactic latitude in degrees
 *
 * @return Galactic latitude in degrees.
 ***************************************************************************/
inline
double GSkyDir::b_deg(void) const
{
    return (b() * gammalib::rad2deg);
}


/***********************************************************************//**
 * @brief Return Right Ascension in radians
 *
 * @return Right Ascension in radians.
 ***************************************************************************/
inline
const double& GSkyDir::ra(void) const
{
    if (!m_has_radec && m_has_lb) {
        gal2equ();
    }
    return m_ra;
}


/***********************************************************************//**
 * @brief Returns Right Ascension in degrees
 *
 * @return Right Ascension in degrees.
 ***************************************************************************/
inline
double GSkyDir::ra_deg(void) const
{
    return (ra() * gammalib::rad2deg);
}


/***********************************************************************//**
 * @brief Return Declination in radians
 *
 * @return Declination in radians.
 ***************************************************************************/
inline
const double& GSkyDir::dec(void) const
{
    if (!m_has_radec && m_has_lb) {
        gal2equ();
    }
    return m_dec;
}


/***********************************************************************//**
 * @brief Returns Declination in degrees
 *
 * @return Declination in degrees.
 ***************************************************************************/
inline
double GSkyDir::dec_deg(void) const
{
    return (dec() * gammalib::rad2deg);
}


/***********************************************************************//**
 * @brief Compute angular distance between sky directions in radians
 *
 * @param[in] dir Sky direction to which distance is to be computed.
 * @return Angular distance (radians).
 *
 * Computes the angular distance between two sky directions in radians.
 ***************************************************************************/
inline
double GSkyDir::dist(const GSkyDir& dir) const
{
    return (gammalib::acos(cos_dist(dir)));
}


/***********************************************************************//**
 * @brief Compute angular distance between sky directions in degrees
 *
 * @param[in] dir Sky direction to which distance is to be computed.
 * @return Angular distance (degrees).
 *
 * Computes the angular distance between two sky directions in degrees.
 ***************************************************************************/
inline
double GSkyDir::dist_deg(const GSkyDir& dir) const
{
    return (dist(dir) * gammalib::rad2deg);
}


/***********************************************************************//**
 * @brief Compute position angle between sky directions in degrees
 *
 * @param[in] dir Sky direction.
 * @return Position angle in degrees.
 *
 * Returns the position angle of a sky direction in degrees.
 *
 * See the posang() method for more information about the computation of the
 * position angle.
 ***************************************************************************/
inline
double GSkyDir::posang_deg(const GSkyDir& dir) const
{
    return (posang(dir) * gammalib::rad2deg);
}

#endif /* GSKYDIR_HPP */
