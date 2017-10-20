/***************************************************************************
 *                  GNdarray.i - N-dimensional array class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016-2017 by Juergen Knoedlseder                         *
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
 * @file GNdarray.i
 * @brief N-dimensional array class interface definition
 * @author Juergen Knoedlseder
 */

%{
/* Put headers and other declarations here that are needed for compilation */
#include "GNdarray.hpp"
%}


/***********************************************************************//**
 * @class GNdarray
 *
 * @brief N-dimensional array class
 *
 * This class implement a n-dimensional double precision floating point
 * array.
 ***************************************************************************/
class GNdarray : public GBase {

public:
    // Constructors and destructors
    GNdarray(void);
    explicit GNdarray(const int& nx);
    GNdarray(const int& nx, const int& ny);
    GNdarray(const int& nx, const int& ny, const int& nz);
    GNdarray(const std::vector<int>& n);
    GNdarray(const GNdarray& array);
    virtual ~GNdarray(void);

    // Operators
    bool      operator==(const GNdarray& array) const;
    bool      operator!=(const GNdarray& array) const;
    GNdarray& operator+=(const GNdarray& array);
    GNdarray& operator-=(const GNdarray& array);
    GNdarray& operator+=(const double& value);
    GNdarray& operator-=(const double& value);
    GNdarray& operator*=(const double& value);
    GNdarray  operator-(void) const;

    // Methods
    void                    clear(void);
    GNdarray*               clone(void) const;
    std::string             classname(void) const;
    int                     dim() const;
    int                     size(void) const;
    const std::vector<int>& shape(void) const;
    const std::vector<int>& strides(void) const;
    void                    shape(const std::vector<int>& shape);
    double*                 data(void);
};


/***********************************************************************//**
 * @brief GNdarray class extension
 ***************************************************************************/
%extend GNdarray {
    double __getitem__(int GTuple[]) {
        // Handle first the special case of a single index. In Python we want
        // to support accessing a n-dimensional array with a single index so
        // that we can iterate over the array
        if (GTuple[0] == 1) {
            if (GTuple[1] >= 0 && GTuple[1] < self->size()) {
                return ((*self)(GTuple[1]));
            }
            else {
                throw GException::out_of_range("__getitem__(int)", "Array "
                      "index", GTuple[1], self->size());
            }
        }
        if (GTuple[0] != self->dim()) {
            throw GException::invalid_value("__getitem__(int)", "Invalid "
                  "access of "+gammalib::str(self->dim())+"-dimensional "
                  "array with "+gammalib::str(GTuple[0])+"-dimensional "
                  "access operator.");
        }
        int index = 0;
        for (size_t i = 0; i < self->shape().size(); ++i) {
            if (GTuple[i+1] < 0 || GTuple[i+1] >= self->shape()[i]) {
                throw GException::out_of_range("__getitem__(int)", "Dimension "+
                      gammalib::str(i)+" index", GTuple[i+1], self->shape()[i]);
            }
            index += self->strides()[i] * GTuple[i+1];
        }
        return ((*self)(index));
    }
    void __setitem__(int GTuple[], const double& value) {
        if (GTuple[0] != self->dim()) {
            throw GException::invalid_value("__setitem__(int)", "Invalid "
                  "access of "+gammalib::str(self->dim())+"-dimensional "
                  "array with "+gammalib::str(GTuple[0])+"-dimensional "
                  "access operator.");
        }
        int index = 0;
        for (size_t i = 0; i < self->shape().size(); ++i) {
            if (GTuple[i+1] < 0 || GTuple[i+1] >= self->shape()[i]) {
                throw GException::out_of_range("__setitem__(int)", "Dimension "+
                      gammalib::str(i)+" index", GTuple[i+1], self->shape()[i]);
            }
            index += self->strides()[i] * GTuple[i+1];
        }
        (*self)(index) = value;
    }
    GNdarray __add__(const GNdarray &a) {
        return (*self) + a;
    }
    GNdarray __add__(const double &a) {
        return (*self) + a;
    }
    GNdarray __sub__(const GNdarray &a) {
        return (*self) - a;
    }
    GNdarray __sub__(const double &a) {
        return (*self) - a;
    }
    GNdarray __mul__(const double &a) {
        return (*self) * a;
    }
    // Python 2.x
    GNdarray __div__(const double &a) {
        return (*self) / a;
    }
    // Python 3.x
    GNdarray __truediv__(const double &a) {
        return (*self) / a;
    }
    // Python 2.x operator/=
    GNdarray __idiv__(const double& value) {
        self->operator/=(value);
        return (*self);
    }
    // Python 3.x operator/=
    GNdarray __itruediv__(const double& value) {
        self->operator/=(value);
        return (*self);
    }
    int __is__(const GNdarray &a) {
            return (*self) == a;
    }
    GNdarray copy() {
        return (*self);
    }
    double min() {
        return min(*self);
    }
    double max() {
        return max(*self);
    }
    double sum() {
        return sum(*self);
    }
    GNdarray acos() {
        return acos(*self);
    }
    GNdarray acosh() {
        return acosh(*self);
    }
    GNdarray asin() {
        return asin(*self);
    }
    GNdarray asinh() {
        return asinh(*self);
    }
    GNdarray atan() {
        return atan(*self);
    }
    GNdarray atanh() {
        return atanh(*self);
    }
    GNdarray cos() {
        return cos(*self);
    }
    GNdarray cosh() {
        return cosh(*self);
    }
    GNdarray exp() {
        return exp(*self);
    }
    GNdarray abs() {
        return abs(*self);
    }
    GNdarray log() {
        return log(*self);
    }
    GNdarray log10() {
        return log10(*self);
    }
    GNdarray sin() {
        return sin(*self);
    }
    GNdarray sinh() {
        return sinh(*self);
    }
    GNdarray sqrt() {
        return sqrt(*self);
    }
    GNdarray tan() {
        return tan(*self);
    }
    GNdarray tanh() {
        return tanh(*self);
    }
};
