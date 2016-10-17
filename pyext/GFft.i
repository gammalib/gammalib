/***************************************************************************
 *                 GFft.i - Fast Fourier transformation class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2016 by Juergen Knoedlseder                              *
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
 * @file GFft.i
 * @brief Fast Fourier transformation class interface definition
 * @author Juergen Knoedlseder
 */

%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFft.hpp"
%}


/***********************************************************************//**
 * @class GFft
 *
 * @brief Fast Fourier Transformation class
 ***************************************************************************/
class GFft : public GBase {

public:
    // Constructors and destructors
    GFft(void);
    explicit GFft(const GNdarray& array);
    GFft(const GFft& fft);
    virtual ~GFft(void);

    // Operators
    GFft& operator+=(const GFft& fft);
    GFft& operator-=(const GFft& fft);
    GFft& operator*=(const GFft& fft);
    GFft& operator/=(const GFft& fft);
    GFft  operator-(void) const;

    // Methods
    void                    clear(void);
    GFft*                   clone(void) const;
    std::string             classname(void) const;
    int                     dim() const;
    int                     size(void) const;
    const std::vector<int>& shape(void) const;
    const std::vector<int>& strides(void) const;
    void                    forward(const GNdarray& array);
    GNdarray                backward(void) const;
};


/***********************************************************************//**
 * @brief GFft class extension
 ***************************************************************************/
%extend GFft {
    std::complex<double> __getitem__(int GTuple[]) {
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
        for (int i = 0; i < self->shape().size(); ++i) {
            if (GTuple[i+1] < 0 || GTuple[i+1] >= self->shape()[i]) {
                throw GException::out_of_range("__getitem__(int)", "Dimension "+
                      gammalib::str(i)+" index", GTuple[i+1], self->shape()[i]);
            }
            index += self->strides()[i] * GTuple[i+1];
        }
        return ((*self)(index));
    }
    void __setitem__(int GTuple[], const std::complex<double>& value) {
        if (GTuple[0] != self->dim()) {
            throw GException::invalid_value("__setitem__(int)", "Invalid "
                  "access of "+gammalib::str(self->dim())+"-dimensional "
                  "array with "+gammalib::str(GTuple[0])+"-dimensional "
                  "access operator.");
        }
        int index = 0;
        for (int i = 0; i < self->shape().size(); ++i) {
            if (GTuple[i+1] < 0 || GTuple[i+1] >= self->shape()[i]) {
                throw GException::out_of_range("__setitem__(int)", "Dimension "+
                      gammalib::str(i)+" index", GTuple[i+1], self->shape()[i]);
            }
            index += self->strides()[i] * GTuple[i+1];
        }
        (*self)(index) = value;
    }
    GFft __add__(const GFft &a) {
        return (*self) + a;
    }
    GFft __sub__(const GFft &a) {
        return (*self) - a;
    }
    GFft __mul__(const GFft &a) {
        return (*self) * a;
    }
    // Python 2.x
    GFft __div__(const GFft &a) {
        return (*self) / a;
    }
    // Python 3.x
    GFft __truediv__(const GFft &a) {
        return (*self) / a;
    }
    // Python 2.x operator/=
    GFft __idiv__(const GFft& a) {
        self->operator/=(a);
        return (*self);
    }
    // Python 3.x operator/=
    GFft __itruediv__(const GFft& a) {
        self->operator/=(a);
        return (*self);
    }
    GFft copy() {
        return (*self);
    }
};
