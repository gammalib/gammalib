/***************************************************************************
 *            GFitsHeader.hpp - FITS header cards container class          *
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
 * @file GFitsHeader.i
 * @brief FITS header cards container class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitsHeader.hpp"
#include "GTools.hpp"
%}

//%include stl.i


/***********************************************************************//**
 * @class GFitsHeader
 *
 * @brief Interface for FITS header class
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
    GFitsHeaderCard& at(const int& cardno);
    GFitsHeaderCard& at(const std::string& keyname);
    std::string      string(const int& cardno) const;
    std::string      string(const std::string& keyname) const;
    double           real(const int& cardno) const;
    double           real(const std::string& keyname) const;
    int              integer(const int& cardno) const;
    int              integer(const std::string& keyname) const;
    int              size(void) const;
    bool             isempty(void) const;
    GFitsHeaderCard& append(const GFitsHeaderCard& card);
    GFitsHeaderCard& insert(const int& cardno, const GFitsHeaderCard& card);
    GFitsHeaderCard& insert(const std::string& keyname, const GFitsHeaderCard& card);
    void             remove(const int& cardno);
    void             remove(const std::string& keyname);
    void             reserve(const int& num);
    void             extend(const GFitsHeader& header);
    bool             contains(const int& cardno) const;
    bool             contains(const std::string& keyname) const;
    void             load(void* vptr);
    void             save(void* vptr) const;
};


/***********************************************************************//**
 * @brief GFitsHeader class extension
 ***************************************************************************/
%extend GFitsHeader {
    GFitsHeaderCard& __getitem__(const int& cardno) {
        if (cardno >= 0 && cardno < self->size()) {
            return (*self)[cardno];
        }
        else {
            throw GException::out_of_range("__getitem__(int)",
                                           "Header card number",
                                           cardno, self->size());
        }
    }
    GFitsHeaderCard& __getitem__(const std::string& keyname) {
        return (*self)[keyname];
    }
    void __setitem__(const int& cardno, const GFitsHeaderCard& card) {
        if (cardno >= 0 && cardno < self->size()) {
            (*self)[cardno] = card;
            return;
        }
        else {
            throw GException::out_of_range("__setitem__(int)",
                                           "Header card number",
                                           cardno, self->size());
        }
    }
    void __setitem__(const std::string& keyname, const GFitsHeaderCard& card) {
        (*self)[keyname] = card;
        return;
    }
    GFitsHeader copy() {
        return (*self);
    }
}
