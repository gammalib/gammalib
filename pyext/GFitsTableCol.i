/***************************************************************************
 *    GFitsTableCol.i  - FITS table column abstract base class SWIG file   *
 * ----------------------------------------------------------------------- *
 *  copyright : (C) 2008-2011 by Jurgen Knodlseder                         *
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
 * @file GFitsTableCol.i
 * @brief FITS table column abstract Python base class definition
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitsTableCol.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GFitsTableCol
 *
 * @brief FITS table column abstract Python base class definition
 ***************************************************************************/
class GFitsTableCol {
public:
    // Constructors and destructors
    GFitsTableCol(void);
    explicit GFitsTableCol(const std::string& name, const int& length,
                           const int& number,       const int& width);
    GFitsTableCol(const GFitsTableCol& column);
    virtual ~GFitsTableCol(void);

    // Virtual Methods
    virtual std::string string(const int& row, const int& inx = 0) = 0;
    virtual double      real(const int& row, const int& inx = 0) = 0;
    virtual int         integer(const int& row, const int& inx = 0) = 0;
    virtual void        insert(const int& rownum, const int& nrows) = 0;
    virtual void        remove(const int& rownum, const int& nrows) = 0;

    // Base class Methods
    void        name(const std::string& name);
    std::string name(void) const;
    int         colnum(void) const;
    int         type(void) const;
    int         repeat(void) const;
    int         width(void) const;
    int         number(void) const;
    int         length(void) const;
    int         anynul(void) const;
};


/***********************************************************************//**
 * @brief GFitsTableCol class extension
 ***************************************************************************/
%extend GFitsTableCol {
    char *__str__() {
        return tochar(self->print());
    }
};
