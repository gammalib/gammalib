/***************************************************************************
 *                  GTimeReference.i - Time reference class                *
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
 * @file GTimeReference.i
 * @brief Time reference class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GTimeReference.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GTimeReference
 *
 * @brief Implements a time reference
 ***************************************************************************/
class GTimeReference : public GBase {
public:
    // Constructors and destructors
    GTimeReference(void);
    GTimeReference(const GTimeReference& ref);
    explicit GTimeReference(const double&      mrdref,
                            const std::string& timeunit,
                            const std::string& timesys = "TT",
                            const std::string& timeref = "local");
    explicit GTimeReference(const int&         mjdrefi,
                            const double&      mrdreff,
                            const std::string& timeunit,
                            const std::string& timesys = "TT",
                            const std::string& timeref = "local");
    explicit GTimeReference(const GFitsHDU* hdu);
    virtual ~GTimeReference(void);

    // Methods
    void               clear(void);
    GTimeReference*    clone(void) const;
    void               read(const GFitsHDU* hdu);
    void               write(GFitsHDU* hdu) const;
    void               set(const double&      mrdref,
                           const std::string& timeunit,
                           const std::string& timesys = "TT",
                           const std::string& timeref = "local");
    void               set(const int&         mjdrefi,
                           const double&      mjdreff,
                           const std::string& timeunit,
                           const std::string& timesys = "TT",
                           const std::string& timeref = "local");
    const double&      mjdref(void) const;
    int                mjdrefi(void) const;
    double             mjdreff(void) const;
    const std::string& timeunit(void) const;
    const std::string& timesys(void) const;
    const std::string& timeref(void) const;
    double             unitseconds(void) const;
};


/***********************************************************************//**
 * @brief GTimeReference class extension
 ***************************************************************************/
%extend GTimeReference {
    GTimeReference copy() {
        return (*self);
    }
};
