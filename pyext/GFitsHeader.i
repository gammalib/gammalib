/***************************************************************************
 *           GFitsHeader.i  - FITS header handling class SWIG file         *
 * ----------------------------------------------------------------------- *
 *  copyright : (C) 2008-2011 by Juergen Knoedlseder                       *
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
 * @file GFitsHeader.i
 * @brief FITS header class Python interface definition
 * @author J. Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitsHeader.hpp"
#include "GTools.hpp"
%}

%include stl.i


/***********************************************************************//**
 * @class GFitsHeader
 *
 * @brief Implements FITS header class SWIG interface
 *
 * The FITS header class contains all cards that are found in the header of
 * a HDU. All cards will be hold in memory, so no link to a FITS file is
 * required. Cards may be read from one file (using the 'open' method) and
 * saved into another file (using the 'save' method). Cards are added or
 * changed using the 'update' method.
 ***************************************************************************/
class GFitsHeader : public GBase {
public:
    // Constructors and destructors
    GFitsHeader(void);
    GFitsHeader(const GFitsHeader& header);
    virtual ~GFitsHeader(void);

    // Methods
    void             clear(void);
    GFitsHeader*     clone(void) const;
    int              size(void) const;
    void             open(void* vptr);
    void             save(void* vptr);
    void             close(void);
    bool             hascard(const std::string& keyname) const;
    bool             hascard(const int& cardno) const;
    void             update(const GFitsHeaderCard& card);
    GFitsHeaderCard* card(const std::string& keyname);
    GFitsHeaderCard* card(const int& cardno);
    std::string      string(const std::string& keyname);
    std::string      string(const int& cardno);
    double           real(const std::string& keyname);
    double           real(const int& cardno);
    int              integer(const std::string& keyname);
    int              integer(const int& cardno);
};


/***********************************************************************//**
 * @brief GFitsHeader class SWIG extension
 ***************************************************************************/
%extend GFitsHeader {
    char *__str__() {
        return tochar(self->print());
    }
    GFitsHeader copy() {
        return (*self);
    }
}
