/***************************************************************************
 *         GFitsTableCol.i  - FITS table column abstract base class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2011 by Jurgen Knodlseder                           *
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
#include "GFitsTableBitCol.hpp"
#include "GFitsTableBoolCol.hpp"
#include "GFitsTableByteCol.hpp"
#include "GFitsTableDoubleCol.hpp"
#include "GFitsTableFloatCol.hpp"
#include "GFitsTableLongCol.hpp"
#include "GFitsTableLongLongCol.hpp"
#include "GFitsTableShortCol.hpp"
#include "GFitsTableStringCol.hpp"
#include "GFitsTableULongCol.hpp"
#include "GFitsTableUShortCol.hpp"
#include "GTools.hpp"
#include "GException.hpp"
%}
%include "std_vector.i"
%template(vectori) std::vector<int>;


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

    // Pure virtual Methods
    virtual void           clear(void) = 0;
    virtual GFitsTableCol* clone(void) const = 0;
    virtual std::string    string(const int& row, const int& inx = 0) const = 0;
    virtual double         real(const int& row, const int& inx = 0) const = 0;
    virtual int            integer(const int& row, const int& inx = 0) const = 0;
    virtual void           insert(const int& rownum, const int& nrows) = 0;
    virtual void           remove(const int& rownum, const int& nrows) = 0;

    // Base class Methods
    void             name(const std::string& name);
    void             unit(const std::string& unit);
    void             dim(const std::vector<int>& dim);
    std::string      name(void) const;
    std::string      unit(void) const;
    std::vector<int> dim(void) const;
    int              colnum(void) const;
    int              type(void) const;
    int              repeat(void) const;
    int              width(void) const;
    int              number(void) const;
    int              length(void) const;
    int              anynul(void) const;
};


/***********************************************************************//**
 * @brief GFitsTableCol class extension
 *
 * A number of __setitem__ methods are implemented that call the respective
 * method of the derived class. This allows transparent setting column
 * elements through the GFitsTableCol class (without any need for type
 * casting). Note that this code produces some warnings of the type
 *
 *    GFitsTableCol.i:115: Warning 467: Overloaded method
 *    GFitsTableCol::__setitem__(int [],std::string) not supported
 *    (no type checking rule for 'int []').
 *
 * on swig 2.0.2, yet it works (as can be verified by checking the
 * wrapper file). A similar problem has been reported on the web here:
 * http://permalink.gmane.org/gmane.comp.programming.swig/16475
 *
 * @todo Add range checking for type casting
 ***************************************************************************/
%extend GFitsTableCol {
    char *__str__() {
        return tochar(self->print());
    }
    
    // String setting
    void __setitem__(int GFitsTableColInx[], std::string value) {
        if (dynamic_cast<GFitsTableStringCol*>(self) != NULL) {
            GFitsTableStringCol* col = dynamic_cast<GFitsTableStringCol*>(self);
            if (GFitsTableColInx[0] == 1)
                (*col)(GFitsTableColInx[1]) = value;
            else
                (*col)(GFitsTableColInx[1], GFitsTableColInx[2]) = value;
        }
        else {
            throw GException::bad_type("GFitsTableCol::__setitem__",
                  "Column type does not support string setting.");            
        }
    }
    
    // Integer setting
    void __setitem__(int GFitsTableColInx[], int value) {
        if (dynamic_cast<GFitsTableBitCol*>(self) != NULL) {
            GFitsTableBitCol* col = dynamic_cast<GFitsTableBitCol*>(self);
            if (GFitsTableColInx[0] == 1)
                (*col)(GFitsTableColInx[1]) = (bool)value;
            else
                (*col)(GFitsTableColInx[1], GFitsTableColInx[2]) = (bool)value;
        }
        else if (dynamic_cast<GFitsTableBoolCol*>(self) != NULL) {
            GFitsTableBoolCol* col = dynamic_cast<GFitsTableBoolCol*>(self);
            if (GFitsTableColInx[0] == 1)
                (*col)(GFitsTableColInx[1]) = (bool)value;
            else
                (*col)(GFitsTableColInx[1], GFitsTableColInx[2]) = (bool)value;
        }
        else if (dynamic_cast<GFitsTableByteCol*>(self) != NULL) {
            GFitsTableByteCol* col = dynamic_cast<GFitsTableByteCol*>(self);
            if (GFitsTableColInx[0] == 1)
                (*col)(GFitsTableColInx[1]) = (unsigned char)value;
            else
                (*col)(GFitsTableColInx[1], GFitsTableColInx[2]) = (unsigned char)value;
        }
        else if (dynamic_cast<GFitsTableLongCol*>(self) != NULL) {
            GFitsTableLongCol* col = dynamic_cast<GFitsTableLongCol*>(self);
            if (GFitsTableColInx[0] == 1)
                (*col)(GFitsTableColInx[1]) = (long)value;
            else
                (*col)(GFitsTableColInx[1], GFitsTableColInx[2]) = (long)value;
        }
        else if (dynamic_cast<GFitsTableLongLongCol*>(self) != NULL) {
            GFitsTableLongLongCol* col = dynamic_cast<GFitsTableLongLongCol*>(self);
            if (GFitsTableColInx[0] == 1)
                (*col)(GFitsTableColInx[1]) = (long long)value;
            else
                (*col)(GFitsTableColInx[1], GFitsTableColInx[2]) = (long long)value;
        }
        else if (dynamic_cast<GFitsTableShortCol*>(self) != NULL) {
            GFitsTableShortCol* col = dynamic_cast<GFitsTableShortCol*>(self);
            if (GFitsTableColInx[0] == 1)
                (*col)(GFitsTableColInx[1]) = (short)value;
            else
                (*col)(GFitsTableColInx[1], GFitsTableColInx[2]) = (short)value;
        }
        else if (dynamic_cast<GFitsTableULongCol*>(self) != NULL) {
            GFitsTableULongCol* col = dynamic_cast<GFitsTableULongCol*>(self);
            if (GFitsTableColInx[0] == 1)
                (*col)(GFitsTableColInx[1]) = (unsigned long)value;
            else
                (*col)(GFitsTableColInx[1], GFitsTableColInx[2]) = (unsigned long)value;
        }
        else if (dynamic_cast<GFitsTableUShortCol*>(self) != NULL) {
            GFitsTableUShortCol* col = dynamic_cast<GFitsTableUShortCol*>(self);
            if (GFitsTableColInx[0] == 1)
                (*col)(GFitsTableColInx[1]) = (unsigned short)value;
            else
                (*col)(GFitsTableColInx[1], GFitsTableColInx[2]) = (unsigned short)value;
        }
        else if (dynamic_cast<GFitsTableDoubleCol*>(self) != NULL) {
            GFitsTableDoubleCol* col = dynamic_cast<GFitsTableDoubleCol*>(self);
            if (GFitsTableColInx[0] == 1)
                (*col)(GFitsTableColInx[1]) = (double)value;
            else
                (*col)(GFitsTableColInx[1], GFitsTableColInx[2]) = (double)value;
        }
        else if (dynamic_cast<GFitsTableFloatCol*>(self) != NULL) {
            GFitsTableFloatCol* col = dynamic_cast<GFitsTableFloatCol*>(self);
            if (GFitsTableColInx[0] == 1)
                (*col)(GFitsTableColInx[1]) = (float)value;
            else
                (*col)(GFitsTableColInx[1], GFitsTableColInx[2]) = (float)value;
        }
        else {
            throw GException::bad_type("GFitsTableCol::__setitem__",
                  "Column type does not support integer setting.");            
        }
    }
    
    // Floating point setting
    void __setitem__(int GFitsTableColInx[], double value) {
        if (dynamic_cast<GFitsTableDoubleCol*>(self) != NULL) {
            GFitsTableDoubleCol* col = dynamic_cast<GFitsTableDoubleCol*>(self);
            if (GFitsTableColInx[0] == 1)
                (*col)(GFitsTableColInx[1]) = value;
            else
                (*col)(GFitsTableColInx[1], GFitsTableColInx[2]) = value;
        }
        else if (dynamic_cast<GFitsTableFloatCol*>(self) != NULL) {
            GFitsTableFloatCol* col = dynamic_cast<GFitsTableFloatCol*>(self);
            if (GFitsTableColInx[0] == 1)
                (*col)(GFitsTableColInx[1]) = (float)value;
            else
                (*col)(GFitsTableColInx[1], GFitsTableColInx[2]) = (float)value;
        }

        else {
            throw GException::bad_type("GFitsTableCol::__setitem__",
                  "Column type does not support floating point setting.");            
        }
    }
    
};
