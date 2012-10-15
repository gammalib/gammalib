/***************************************************************************
 *          GNodeArray.i  -  Array of nodes class SWIG definition          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2012 by Juergen Knoedlseder                         *
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
 * @file GNodeArray.i
 * @brief Node array class Python interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GNodeArray.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GNodeArray
 *
 * @brief Python interface for the node array class
 ***************************************************************************/
class GNodeArray : public GBase {

public:
    // Constructors and destructors
    GNodeArray(void);
    GNodeArray(const GNodeArray& array);
    GNodeArray(const int& num, const double* array);
    GNodeArray(const GVector& vector);
    GNodeArray(const std::vector<double>& vector);
    virtual ~GNodeArray(void);

    // Methods
    void          clear(void);
    GNodeArray*   clone(void) const;
    int           size(void) const;
    void          nodes(const int& num, const double* array);
    void          nodes(const GVector& vector);
    void          nodes(const std::vector<double>& vector);
    void          append(const double& node);
    double        interpolate(const double& value,
                              const std::vector<double>& vector) const;
    void          set_value(const double& value);
    const int&    inx_left(void) const;
    const int&    inx_right(void) const;
    const double& wgt_left(void) const;
    const double& wgt_right(void) const;
};


/***********************************************************************//**
 * @brief GNodeArray class extension
 ***************************************************************************/
%extend GNodeArray {
    char *__str__() {
        return tochar(self->print());
    }
    double __getitem__(const int& index) {
        if (index >= 0 && index < self->size())
            return (*self)[index];
        else
            throw GException::out_of_range("__getitem__(int)", index, self->size());
    }
    void __setitem__(const int& index, const double& val) {
        if (index>=0 && index < self->size()) {
            (*self)[index] = val;
            return;
        }
        else
            throw GException::out_of_range("__setitem__(int)", index, self->size());
    }
    GNodeArray copy() {
        return (*self);
    }
};
