/***************************************************************************
 *                       GVector.i - Vector class                          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2018 by Juergen Knoedlseder                         *
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
 * @file GVector.i
 * @brief Vector class Python interface definition
 * @author Juergen Knoedlseder
 */

%{
/* Put headers and other declarations here that are needed for compilation */
#include "GVector.hpp"
#include "GException.hpp"
%}


/***********************************************************************//**
 * @class GVector
 *
 * @brief Vector class
 ***************************************************************************/
class GVector : public GBase {
public:
    // Constructors and destructors
    GVector(void);
    explicit GVector(const int& num);
    explicit GVector(const double& a);
    explicit GVector(const double& a, const double& b);
    explicit GVector(const double& a, const double& b, const double& c);
    GVector(const GVector& v);
    virtual ~GVector(void);

    // Vector operators
    bool     operator==(const GVector& vector) const;
    bool     operator!=(const GVector& vector) const;
    GVector& operator+=(const GVector& vector);
    GVector& operator-=(const GVector& vector);
    GVector& operator+=(const double& value);
    GVector& operator-=(const double& value);
    GVector& operator*=(const double& value);
    GVector  operator-(void) const;

    // Vector functions
    void        clear(void);
    GVector*    clone(void) const;
    std::string classname(void) const;
    const int&  size(void) const;
    int         non_zeros(void) const;
};


/***********************************************************************//**
 * @brief GVector class extension
 ***************************************************************************/
%extend GVector {
    double __getitem__(const int& index) {
        if (index >= 0 && index < self->size()) {
            return (*self)[index];
        }
        else {
            throw GException::out_of_range("__getitem__(int)", index, self->size());
        }
    }
    void __setitem__(const int& index, const double& val) {
        if (index>=0 && index < self->size()) {
            (*self)[index] = val;
        }
        else {
            throw GException::out_of_range("__setitem__(int)", index, self->size());
        }
    }
    GVector __add__(const GVector &a) {
        return (*self) + a;
    }
    GVector __add__(const double &a) {
        return (*self) + a;
    }
    GVector __sub__(const GVector &a) {
        return (*self) - a;
    }
    GVector __sub__(const double &a) {
        return (*self) - a;
    }
    double __mul__(const GVector &a) {
        return (*self) * a;
    }
    GVector __mul__(const double &a) {
        return (*self) * a;
    }
    // Python 2.x
    GVector __div__(const double &a) {
        return (*self) / a;
    }
    // Python 3.x
    GVector __truediv__(const double &a) {
        return (*self) / a;
    }
    // Python 2.x operator/=
    GVector __idiv__(const double& value) {
        self->operator/=(value);
        return (*self);
    }
    // Python 3.x operator/=
    GVector __itruediv__(const double& value) {
        self->operator/=(value);
        return (*self);
    }
    int __is__(const GVector &a) {
            return (*self) == a;
    }
    GVector copy() {
        return (*self);
    }
    GVector cross(const GVector &a) {
        return cross(*self, a);
    }
    double norm() {
        return norm(*self);
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
    GVector acos() {
        return acos(*self);
    }
    GVector acosh() {
        return acosh(*self);
    }
    GVector asin() {
        return asin(*self);
    }
    GVector asinh() {
        return asinh(*self);
    }
    GVector atan() {
        return atan(*self);
    }
    GVector atanh() {
        return atanh(*self);
    }
    GVector cos() {
        return cos(*self);
    }
    GVector cosh() {
        return cosh(*self);
    }
    GVector exp() {
        return exp(*self);
    }
    GVector abs() {
        return abs(*self);
    }
    GVector log() {
        return log(*self);
    }
    GVector log10() {
        return log10(*self);
    }
    GVector sin() {
        return sin(*self);
    }
    GVector sinh() {
        return sinh(*self);
    }
    GVector sqrt() {
        return sqrt(*self);
    }
    GVector tan() {
        return tan(*self);
    }
    GVector tanh() {
        return tanh(*self);
    }
};
