/***************************************************************************
 *             GCTAInstDir.hpp - CTA instrument direction class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2020 by Juergen Knoedlseder                         *
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
#include <cstdint>
#include "GInstDir.hpp"
#include "GSkyDir.hpp"
#include "GMath.hpp"
#include "GException.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_DETX                                          "GCTAInstDir::detx()"
#define G_DETY                                          "GCTAInstDir::dety()"
#define G_THETA                                        "GCTAInstDir::theta()"
#define G_PHI                                          "GCTAInstDir::theta()"
#define G_DIR                                            "GCTAInstDir::dir()"


/***********************************************************************//**
 * @class GCTAInstDir
 *
 * @brief CTA instrument direction class.
 *
 * The CTA instrument direction defines the measured (or reconstructed)
 * direction of an event. The instrument direction comprises a sky direction
 * and detector coordinates. Not all the information is necessary provided,
 * and the has_dir(), has_detx() and has_dety() methods inform which
 * information is contained in a class instance. Note also that the dir(),
 * detx() and dety() methods will throw an exception in case that no
 * information exists. This assures that only valid information can be
 * accessed.
 ***************************************************************************/
class GCTAInstDir : public GInstDir {

public:
    // Constructors and destructors
    GCTAInstDir(void);
    explicit GCTAInstDir(const GSkyDir& dir);
    GCTAInstDir(const double& detx, const double& dety);
    GCTAInstDir(const GSkyDir& dir, const double& detx, const double& dety);
    GCTAInstDir(const GCTAInstDir& dir);
    virtual ~GCTAInstDir(void);

    // Operators
    GCTAInstDir& operator=(const GCTAInstDir& dir);

    // Implemented pure virtual base class methods
    virtual void         clear(void);
    virtual GCTAInstDir* clone(void) const;
    virtual std::string  classname(void) const;
    virtual uint64_t     hash(void) const;
    virtual std::string  print(const GChatter& chatter = NORMAL) const;

    // Other methods
    void           dir(const GSkyDir& dir);
    const GSkyDir& dir(void) const;
    void           detx(const double &x);
    void           dety(const double &y);
    const double&  detx(void) const;
    const double&  dety(void) const;
    double         theta(void) const;
    double         phi(void) const;
    const bool&    has_dir(void) const;
    const bool&    has_detx(void) const;
    const bool&    has_dety(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCTAInstDir& dir);
    void free_members(void);

    // Data members
    GSkyDir m_dir;      //!< Observed incident direction of event
    double  m_detx;     //!< Instrument coordinate X (radians)
    double  m_dety;     //!< Instrument coordinate Y (radians)
    bool    m_has_dir;  //!< Has valid incident direction
    bool    m_has_detx; //!< Has valid instrument coordinate X
    bool    m_has_dety; //!< Has valid instrument coordinate Y
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCTAInstDir").
 ***************************************************************************/
inline
std::string GCTAInstDir::classname(void) const
{
    return ("GCTAInstDir");
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
    m_dir     = dir;
    m_has_dir = true;
    return;
}


/***********************************************************************//**
 * @brief Return reference to sky direction (const version)
 *
 * @return Reference to sky direction.
 *
 * @exception GException::runtime_error
 *            Instrument coordinate has no valid sky direction
 *
 * Returns reference to the sky direction.
 ***************************************************************************/
inline
const GSkyDir& GCTAInstDir::dir(void) const
{
    if (!m_has_dir) {
        std::string msg = "Instrument coordinate has no valid sky direction.";
        throw GException::runtime_error(G_DIR, msg);
    }
    return (m_dir);
}


/***********************************************************************//**
 * @brief Set DETX coordinate (in radians)
 *
 * @param[in] x DETX coordinate (in radians).
 *
 * Set DETX coordinate.
 ***************************************************************************/
inline
void GCTAInstDir::detx(const double &x)
{
    m_detx     = x;
    m_has_detx = true;
    return;
}


/***********************************************************************//**
 * @brief Set DETY coordinate (in radians)
 *
 * @param[in] y DETY coordinate (in radians).
 *
 * Set DETY coordinate.
 ***************************************************************************/
inline
void GCTAInstDir::dety(const double &y)
{
    m_dety     = y;
    m_has_dety = true;
    return;
}


/***********************************************************************//**
 * @brief Return reference to DETX coordinate (in radians)
 *
 * @return Reference to DETX coordinate (in radians).
 *
 * @exception GException::runtime_error
 *            Instrument coordinate has no valid DETX coordinate
 *
 * Returns reference DETX coordinate (in radians).
 ***************************************************************************/
inline
const double& GCTAInstDir::detx(void) const
{
    if (!m_has_detx) {
        std::string msg = "Instrument coordinate has no valid DETX coordinate.";
        throw GException::runtime_error(G_DETX, msg);
    }
    return (m_detx);
}


/***********************************************************************//**
 * @brief Return reference to DETY coordinate (in radians)
 *
 * @return Reference to DETY coordinate (in radians).
 *
 * @exception GException::runtime_error
 *            Instrument coordinate has no valid DETY coordinate
 *
 * Returns reference DETY coordinate (in radians).
 ***************************************************************************/
inline
const double& GCTAInstDir::dety(void) const
{
    if (!m_has_dety) {
        std::string msg = "Instrument coordinate has no valid DETY coordinate.";
        throw GException::runtime_error(G_DETY, msg);
    }
    return (m_dety);
}


/***********************************************************************//**
 * @brief Return offset angle (in radians)
 *
 * @return Offset angle (in radians).
 *
 * @exception GException::runtime_error
 *            Instrument coordinate has no valid DETX or DETY coordinate
 *
 * Returns the offset angle from the camera centre (in radians). The offset
 * angle \f$\theta\f$ is computed using
 *
 * \f[\theta = \sqrt{{\rm DETX}^2 + {\rm DETY}^2}\f]
 *
 * where \f${\rm DETX}\f$ and \f${\rm DETY}\f$ are the detector coordinates.
 ***************************************************************************/
inline
double GCTAInstDir::theta(void) const
{
    if (!m_has_detx || !m_has_dety) {
        std::string msg = "Instrument coordinate has no valid DETX or DETY coordinate.";
        throw GException::runtime_error(G_THETA, msg);
    }
    return (std::sqrt(m_detx*m_detx + m_dety*m_dety));
}


/***********************************************************************//**
 * @brief Return azimuth angle (in radians)
 *
 * @return Azimuth angle (in radians).
 *
 * @exception GException::runtime_error
 *            Instrument coordinate has no valid DETX or DETY coordinate
 *
 * Returns the azimuth angle from the camera centre (in radians). The
 * azimuth angle \f$\phi\f$ is computed using
 *
 * \f[\phi = \arctan \left( \frac{{\rm DETY}}{{\rm DETX}} \right)\f]
 *
 * where \f${\rm DETX}\f$ and \f${\rm DETY}\f$ are the detector coordinates.
 * The \c atan2(DETY,DETX) function is used for the computation.
 ***************************************************************************/
inline
double GCTAInstDir::phi(void) const
{
    if (!m_has_detx || !m_has_dety) {
        std::string msg = "Instrument coordinate has no valid DETX or DETY coordinate.";
        throw GException::runtime_error(G_PHI, msg);
    }
    return (gammalib::atan2(m_dety, m_detx));
}


/***********************************************************************//**
 * @brief Signal if instrument direction has valid sky direction
 *
 * @return True of instrument direction has valid sky direction.
 ***************************************************************************/
inline
const bool& GCTAInstDir::has_dir(void) const
{
    return (m_has_dir);
}


/***********************************************************************//**
 * @brief Signal if instrument direction has valid DETX coordinate
 *
 * @return True of instrument direction has valid DETX coordinate.
 ***************************************************************************/
inline
const bool& GCTAInstDir::has_detx(void) const
{
    return (m_has_detx);
}


/***********************************************************************//**
 * @brief Signal if instrument direction has valid DETY coordinate
 *
 * @return True of instrument direction has valid DETY coordinate.
 ***************************************************************************/
inline
const bool& GCTAInstDir::has_dety(void) const
{
    return (m_has_dety);
}

#endif /* GCTAINSTDIR_HPP */
