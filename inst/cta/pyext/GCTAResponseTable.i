/***************************************************************************
 *             GCTAResponseTable.i  -  CTA response table class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012 by Juergen Knoedlseder                              *
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
 * @file GCTAResponseTable.i
 * @brief CTA response table class definition
 * @author J. Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GTools.hpp"
#include "GCTAResponseTable.hpp"
%}


/***********************************************************************//**
 * @class GCTAResponse
 *
 * @brief Interface for the CTA response table class
 ***************************************************************************/
class GCTAResponseTable {

public:
    // Constructors and destructors
    GCTAResponseTable(void);
    GCTAResponseTable(const GCTAResponseTable& table);
    GCTAResponseTable(const GFitsTable* hdu);
    virtual ~GCTAResponseTable(void);

    // Methods
    void               clear(void);
    GCTAResponseTable* clone(void) const;
    int                size(void) const;
    int                axis(const int& index) const;
    double             axis_lo(const int& index, const int& bin) const;
    double             axis_hi(const int& index, const int& bin) const;
    void               axis_linear(const int& index);
    void               axis_log10(const int& index);
    void               read(const GFitsTable* hdu);
    void               write(GFitsTable* hdu) const;
};


/***********************************************************************//**
 * @brief GCTAResponse class extension
 ***************************************************************************/
%extend GCTAResponseTable {
    char *__str__() {
        return tochar(self->print());
    }
};
