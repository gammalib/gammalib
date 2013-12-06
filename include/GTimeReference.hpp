/***************************************************************************
 *                 GTimeReference.hpp - Time reference class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2013 by Juergen Knoedlseder                         *
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
 * @file GTimeReference.hpp
 * @brief Time reference class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GTIMEREFERENCE_HPP
#define GTIMEREFERENCE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GFitsHDU.hpp"


/***********************************************************************//**
 * @class GTimeReference
 *
 * @brief Implements a time reference
 *
 * This class implements the reference of a time with respect to Modified
 * Julian Days. The zero point of the time is specified by the member
 * m_mjdref that is given in Modified Julian Days. A time may be either
 * either given in seconds or in days. Furthermore, the time system needs to
 * be specified. So far, only Terrestrial Time (TT) is supported. Also the
 * reference of the time system needs to be given. So far, only Local is
 * supported.
 *
 * References:
 * http://aa.usno.navy.mil/publications/docs/Circular_179.php
 * http://hpiers.obspm.fr/eop-pc/index.php?index=TAI-UTC_tab&lang=en
 ***************************************************************************/
class GTimeReference : public GBase {

public:
    // Constructors and destructors
    GTimeReference(void);
    GTimeReference(const GTimeReference& ref);
    explicit GTimeReference(const double&      mrdref,
                            const std::string& timeunit,
                            const std::string& timesys = "TT",
                            const std::string& timeref = "LOCAL");
    explicit GTimeReference(const int&         mjdrefi,
                            const double&      mrdreff,
                            const std::string& timeunit,
                            const std::string& timesys = "TT",
                            const std::string& timeref = "LOCAL");
    explicit GTimeReference(const GFitsHDU& hdu);
    virtual ~GTimeReference(void);
 
    // Operators
    GTimeReference& operator=(const GTimeReference& ref);

    // Methods
    void               clear(void);
    GTimeReference*    clone(void) const;
    void               read(const GFitsHDU& hdu);
    void               write(GFitsHDU& hdu) const;
    void               set(const double&      mrdref,
                           const std::string& timeunit,
                           const std::string& timesys = "TT",
                           const std::string& timeref = "LOCAL");
    void               set(const int&         mjdrefi,
                           const double&      mjdreff,
                           const std::string& timeunit,
                           const std::string& timesys = "TT",
                           const std::string& timeref = "LOCAL");
    const double&      mjdref(void) const;
    int                mjdrefi(void) const;
    double             mjdreff(void) const;
    const std::string& timeunit(void) const;
    const std::string& timesys(void) const;
    const std::string& timeref(void) const;
    double             unitseconds(void) const;
    std::string        print(const GChatter& chatter = NORMAL) const;
  
protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GTimeReference& ref);
    void free_members(void);

    // Protected data members
    double      m_mjdref;        //!< Time MJD reference (days)
    std::string m_timeunit;      //!< Time unit
    std::string m_timesys;       //!< Time system
    std::string m_timeref;       //!< Time reference
    bool        m_unit_sec;      //!< True: unit is seconds, False: unit is days
};

#endif /* GTIMEREFERENCE_HPP */
