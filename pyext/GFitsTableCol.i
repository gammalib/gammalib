/***************************************************************************
 *         GFitsTableCol.i - FITS abstract table column base class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2021 by Juergen Knoedlseder                         *
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
 * @brief FITS abstract table column base class definition
 * @author Juergen Knoedlseder
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


/***********************************************************************//**
 * @class GFitsTableCol
 *
 * @brief FITS table column abstract Python base class definition
 ***************************************************************************/
class GFitsTableCol : public GBase {

public:
    // Constructors and destructors
    GFitsTableCol(void);
    GFitsTableCol(const std::string& name,
                  const int&         nrows,
                  const int&         number,
                  const int&         width);
    GFitsTableCol(const GFitsTableCol& column);
    virtual ~GFitsTableCol(void);

    // Pure virtual methods
    virtual void            clear(void) = 0;
    virtual GFitsTableCol*  clone(void) const = 0;
    virtual std::string     classname(void) const = 0;
    virtual std::string     string(const int& row, const int& inx = 0) const = 0;
    virtual double          real(const int& row, const int& inx = 0) const = 0;
    virtual int             integer(const int& row, const int& inx = 0) const = 0;
    virtual void            insert(const int& row, const int& nrows) = 0;
    virtual void            remove(const int& row, const int& nrows) = 0;
    virtual bool            is_loaded(void) const = 0;

    // Other methods
    void                    name(const std::string& name);
    const std::string&      name(void) const;
    void                    unit(const std::string& unit);
    const std::string&      unit(void) const;
    void                    dim(const std::vector<int>& dim);
    const std::vector<int>& dim(void) const;
    void                    colnum(const int& colnum);
    const int&              colnum(void) const;
    void                    type(const int& type);
    const int&              type(void) const;
    void                    repeat(const int& repeat);
    const int&              repeat(void) const;
    void                    width(const int& width);
    const int&              width(void) const;
    void                    number(const int& number);
    const int&              number(void) const;
    void                    elements(const int& row, const int& elements);
    int                     elements(const int& row) const;
    void                    nrows(const int& nrows);
    const int&              nrows(void) const;
    void                    is_variable(const bool& variable);
    const bool&             is_variable(void) const;
    void                    anynul(const int& anynul);
    const int&              anynul(void) const;
    void                    tscale(const double& tscale);
    const double&           tscale(void) const;
    std::string             tform_binary(void) const;
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

    // String setting
    void __setitem__(int GTuple1D2D[], std::string value) {
        if (dynamic_cast<GFitsTableStringCol*>(self) != NULL) {
            GFitsTableStringCol* col = dynamic_cast<GFitsTableStringCol*>(self);
            if (GTuple1D2D[0] == 1)
                (*col)(GTuple1D2D[1]) = value;
            else
                (*col)(GTuple1D2D[1], GTuple1D2D[2]) = value;
        }
        else {
            throw GException::invalid_argument("GFitsTableCol::__setitem__",
                  "Column type does not support string setting.");
        }
    }

    // Integer setting
    void __setitem__(int GTuple1D2D[], int value) {
        if (dynamic_cast<GFitsTableBitCol*>(self) != NULL) {
            GFitsTableBitCol* col = dynamic_cast<GFitsTableBitCol*>(self);
            if (GTuple1D2D[0] == 1)
                (*col)(GTuple1D2D[1]) = (bool)value;
            else
                (*col)(GTuple1D2D[1], GTuple1D2D[2]) = (bool)value;
        }
        else if (dynamic_cast<GFitsTableBoolCol*>(self) != NULL) {
            GFitsTableBoolCol* col = dynamic_cast<GFitsTableBoolCol*>(self);
            if (GTuple1D2D[0] == 1)
                (*col)(GTuple1D2D[1]) = (bool)value;
            else
                (*col)(GTuple1D2D[1], GTuple1D2D[2]) = (bool)value;
        }
        else if (dynamic_cast<GFitsTableByteCol*>(self) != NULL) {
            GFitsTableByteCol* col = dynamic_cast<GFitsTableByteCol*>(self);
            if (GTuple1D2D[0] == 1)
                (*col)(GTuple1D2D[1]) = (unsigned char)value;
            else
                (*col)(GTuple1D2D[1], GTuple1D2D[2]) = (unsigned char)value;
        }
        else if (dynamic_cast<GFitsTableLongCol*>(self) != NULL) {
            GFitsTableLongCol* col = dynamic_cast<GFitsTableLongCol*>(self);
            if (GTuple1D2D[0] == 1)
                (*col)(GTuple1D2D[1]) = (long)value;
            else
                (*col)(GTuple1D2D[1], GTuple1D2D[2]) = (long)value;
        }
        else if (dynamic_cast<GFitsTableLongLongCol*>(self) != NULL) {
            GFitsTableLongLongCol* col = dynamic_cast<GFitsTableLongLongCol*>(self);
            if (GTuple1D2D[0] == 1)
                (*col)(GTuple1D2D[1]) = (long long)value;
            else
                (*col)(GTuple1D2D[1], GTuple1D2D[2]) = (long long)value;
        }
        else if (dynamic_cast<GFitsTableShortCol*>(self) != NULL) {
            GFitsTableShortCol* col = dynamic_cast<GFitsTableShortCol*>(self);
            if (GTuple1D2D[0] == 1)
                (*col)(GTuple1D2D[1]) = (short)value;
            else
                (*col)(GTuple1D2D[1], GTuple1D2D[2]) = (short)value;
        }
        else if (dynamic_cast<GFitsTableULongCol*>(self) != NULL) {
            GFitsTableULongCol* col = dynamic_cast<GFitsTableULongCol*>(self);
            if (GTuple1D2D[0] == 1)
                (*col)(GTuple1D2D[1]) = (unsigned long)value;
            else
                (*col)(GTuple1D2D[1], GTuple1D2D[2]) = (unsigned long)value;
        }
        else if (dynamic_cast<GFitsTableUShortCol*>(self) != NULL) {
            GFitsTableUShortCol* col = dynamic_cast<GFitsTableUShortCol*>(self);
            if (GTuple1D2D[0] == 1)
                (*col)(GTuple1D2D[1]) = (unsigned short)value;
            else
                (*col)(GTuple1D2D[1], GTuple1D2D[2]) = (unsigned short)value;
        }
        else if (dynamic_cast<GFitsTableDoubleCol*>(self) != NULL) {
            GFitsTableDoubleCol* col = dynamic_cast<GFitsTableDoubleCol*>(self);
            if (GTuple1D2D[0] == 1)
                (*col)(GTuple1D2D[1]) = (double)value;
            else
                (*col)(GTuple1D2D[1], GTuple1D2D[2]) = (double)value;
        }
        else if (dynamic_cast<GFitsTableFloatCol*>(self) != NULL) {
            GFitsTableFloatCol* col = dynamic_cast<GFitsTableFloatCol*>(self);
            if (GTuple1D2D[0] == 1)
                (*col)(GTuple1D2D[1]) = (float)value;
            else
                (*col)(GTuple1D2D[1], GTuple1D2D[2]) = (float)value;
        }
        else {
            throw GException::invalid_argument("GFitsTableCol::__setitem__",
                  "Column type does not support integer setting.");
        }
    }

    // Floating point setting
    void __setitem__(int GTuple1D2D[], double value) {
        if (dynamic_cast<GFitsTableDoubleCol*>(self) != NULL) {
            GFitsTableDoubleCol* col = dynamic_cast<GFitsTableDoubleCol*>(self);
            if (GTuple1D2D[0] == 1)
                (*col)(GTuple1D2D[1]) = value;
            else
                (*col)(GTuple1D2D[1], GTuple1D2D[2]) = value;
        }
        else if (dynamic_cast<GFitsTableFloatCol*>(self) != NULL) {
            GFitsTableFloatCol* col = dynamic_cast<GFitsTableFloatCol*>(self);
            if (GTuple1D2D[0] == 1)
                (*col)(GTuple1D2D[1]) = (float)value;
            else
                (*col)(GTuple1D2D[1], GTuple1D2D[2]) = (float)value;
        }

        else {
            throw GException::invalid_argument("GFitsTableCol::__setitem__",
                  "Column type does not support floating point setting.");
        }
    }
};
