/***************************************************************************
 *         GFitsImageDouble.i - Double precision FITS image class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2021 by Juergen Knoedlseder                         *
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
 * @file GFitsImageDouble.i
 * @brief Double precision FITS image class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitsImageDouble.hpp"
%}

/* __ Typemaps ___________________________________________________________ */
%typemap(in) const double* pixels {
    if (!PyList_Check($input)) {
        SWIG_exception(SWIG_ValueError, "List expected");
    }
    Py_ssize_t num = PyList_Size($input);
    if (num > 0) {
        $1 = (double*)malloc(num*sizeof(double));
        for (int i = 0; i < num; ++i) {
            PyObject *s = PyList_GetItem($input,i);
            if (!PyFloat_Check(s)) {
                free($1);
                SWIG_exception(SWIG_ValueError, "List items must be floats");
            }
            $1[i] = (double)PyFloat_AsDouble(s);
        }
    }
}
%typemap(typecheck, precedence=SWIG_TYPECHECK_INTEGER) const double* pixels {
    $1 = PyList_Check($input) ? 1 : 0;
}
%typemap(freearg) const double* pixels {
    if ($1) {
        free($1);
    }
}


/***********************************************************************//**
 * @class GFitsImageDouble
 *
 * @brief Double precision FITS image class
 ***************************************************************************/
class GFitsImageDouble : public GFitsImage {
public:
    // Constructors and destructors
    GFitsImageDouble(void);
    GFitsImageDouble(const int& nx, const double* pixels = NULL);
    GFitsImageDouble(const int& nx, const int& ny, const double* pixels = NULL);
    GFitsImageDouble(const int& nx, const int& ny, const int& nz, const double* pixels = NULL);
    GFitsImageDouble(const int& nx, const int& ny, const int& nz, const int& nt, const double* pixels = NULL);
    GFitsImageDouble(const std::vector<int>& naxes, const double* pixels = NULL);
    GFitsImageDouble(const GFitsImage& image);
    GFitsImageDouble(const GFitsImageDouble& image);
    virtual ~GFitsImageDouble(void);

    // Methods
    void              clear(void);
    GFitsImageDouble* clone(void) const;
    std::string       classname(void) const;
    double            pixel(const int& ix) const;
    double            pixel(const int& ix, const int& iy) const;
    double            pixel(const int& ix, const int& iy, const int& iz) const;
    double            pixel(const int& ix, const int& iy, const int& iz, const int& it) const;
    int               type(void) const;
};


/***********************************************************************//**
 * @brief GFitsImageDouble class extension
 ***************************************************************************/
%extend GFitsImageDouble {
    double __getitem__(int GTuple1D4D[]) {
        // Check image dimensions
        for (int i = 0; i < GTuple1D4D[0]; ++i) {
             if (GTuple1D4D[i+1] < 0 || GTuple1D4D[i+1] >= self->naxes(i)) {
                throw GException::out_of_range("__getitem__(int)",
                                               "FITS image axis "+
                                               gammalib::str(i)+" index",
                                               GTuple1D4D[i+1],
                                               self->naxes(i));
            }
        }
        // Return pixel value
        if (GTuple1D4D[0] == 1) {
            return (*self)(GTuple1D4D[1]);
        }
        else if (GTuple1D4D[0] == 2) {
            return (*self)(GTuple1D4D[1], GTuple1D4D[2]);
        }
        else if (GTuple1D4D[0] == 3) {
            return (*self)(GTuple1D4D[1], GTuple1D4D[2], GTuple1D4D[3]);
        }
        else if (GTuple1D4D[0] == 4) {
            return (*self)(GTuple1D4D[1], GTuple1D4D[2], GTuple1D4D[3],
                           GTuple1D4D[4]);
        }
        else {
            std::string msg = "Image dimension "+gammalib::str(GTuple1D4D[0])+
                              " must be within [1,4].";
            throw GException::invalid_argument("__getitem__(int)", msg);
        }
    }
    void __setitem__(int GTuple1D4D[], double value) {
        if (GTuple1D4D[0] == 1) {
            self->at(GTuple1D4D[1]) = value;
        }
        else if (GTuple1D4D[0] == 2) {
            self->at(GTuple1D4D[1], GTuple1D4D[2]) = value;
        }
        else if (GTuple1D4D[0] == 3) {
            self->at(GTuple1D4D[1], GTuple1D4D[2], GTuple1D4D[3]) =
                     value;
        }
        else if (GTuple1D4D[0] == 4) {
            self->at(GTuple1D4D[1], GTuple1D4D[2], GTuple1D4D[3],
                     GTuple1D4D[4]) = value;
        }
        else {
            std::string msg = "Image dimension "+gammalib::str(GTuple1D4D[0])+
                              " must be within [1,4].";
            throw GException::invalid_argument("__setitem__(int)", msg);
        }
    }
    PyObject* pixels(void) {
        PyObject *list = PyList_New(self->npix());
        for (int i = 0; i < self->npix(); ++i) {
            PyList_SetItem(list, i, PyFloat_FromDouble((static_cast<double*>(self->pixels())[i])));
        }
        return list;
    }
    GFitsImageDouble copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        state = (gammalib.GFitsHDU.__getstate__(self),
                 [self.naxes(axis) for axis in range(self.naxis())],
                 self.pixels())
        return state
    def __setstate__(self, state):
        self.__init__(state[1], state[2])
        gammalib.GFitsHDU.__setstate__(self, state[0])
}
};
