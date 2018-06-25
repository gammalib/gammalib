/***************************************************************************
 *          GFitsImageFloat.i - Single precision FITS image class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2018 by Juergen Knoedlseder                         *
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
 * @file GFitsImageFloat.i
 * @brief Single precision FITS image class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitsImageFloat.hpp"
%}

/* __ Typemaps ___________________________________________________________ */
%typemap(in) const float* pixels {
    if (!PyList_Check($input)) {
        SWIG_exception(SWIG_ValueError, "List expected");
    }
    int num = PyList_Size($input);
    if (num > 0) {
        $1 = (float*)malloc(num*sizeof(float));
        for (int i = 0; i < num; ++i) {
        PyObject *s = PyList_GetItem($input,i);
        if (!PyFloat_Check(s)) {
            free($1);
            SWIG_exception(SWIG_ValueError, "List items must be floats");
        }
        $1[i] = (float)PyFloat_AsDouble(s);
    }
    }
}
%typemap(typecheck, precedence=SWIG_TYPECHECK_INTEGER) const float* pixels {
    $1 = PyList_Check($input) ? 1 : 0;
}
%typemap(freearg) const float* pixels {
    if ($1) {
        free($1);
    }
}


/***********************************************************************//**
 * @class GFitsImageFloat
 *
 * @brief Single precision FITS image class
 ***************************************************************************/
class GFitsImageFloat : public GFitsImage {
public:
    // Constructors and destructors
    GFitsImageFloat(void);
    GFitsImageFloat(const int& nx, const float* pixels = NULL);
    GFitsImageFloat(const int& nx, const int& ny, const float* pixels = NULL);
    GFitsImageFloat(const int& nx, const int& ny, const int& nz, const float* pixels = NULL);
    GFitsImageFloat(const int& nx, const int& ny, const int& nz, const int& nt, const float* pixels = NULL);
    GFitsImageFloat(const std::vector<int>& naxes, const float* pixels = NULL);
    GFitsImageFloat(const GFitsImage& image);
    GFitsImageFloat(const GFitsImageFloat& image);
    virtual ~GFitsImageFloat(void);

    // Methods
    void             clear(void);
    GFitsImageFloat* clone(void) const;
    std::string      classname(void) const;
    double           pixel(const int& ix) const;
    double           pixel(const int& ix, const int& iy) const;
    double           pixel(const int& ix, const int& iy, const int& iz) const;
    double           pixel(const int& ix, const int& iy, const int& iz, const int& it) const;
    int              type(void) const;
};


/***********************************************************************//**
 * @brief GFitsImageFloat class extension
 ***************************************************************************/
%extend GFitsImageFloat {
    float __getitem__(int GFitsImageInx[]) {
        // Check image dimensions
        for (int i = 0; i < GFitsImageInx[0]; ++i) {
             if (GFitsImageInx[i+1] < 0 || GFitsImageInx[i+1] >= self->naxes(i)) {
                throw GException::out_of_range("__getitem__(int)",
                                               "FITS image axis "+
                                               gammalib::str(i)+" index",
                                               GFitsImageInx[i+1],
                                               self->naxes(i));
            }
        }
        // Return pixel value
        if (GFitsImageInx[0] == 1) {
            return (*self)(GFitsImageInx[1]);
        }
        else if (GFitsImageInx[0] == 2) {
            return (*self)(GFitsImageInx[1], GFitsImageInx[2]);
        }
        else if (GFitsImageInx[0] == 3) {
            return (*self)(GFitsImageInx[1], GFitsImageInx[2], GFitsImageInx[3]);
        }
        else if (GFitsImageInx[0] == 4) {
            return (*self)(GFitsImageInx[1], GFitsImageInx[2], GFitsImageInx[3],
                           GFitsImageInx[4]);
        }
        else {
            throw GException::fits_wrong_image_operator("__getitem__(int)",
                                                        self->naxis(), GFitsImageInx[0]);
        }
    }
    void __setitem__(int GFitsImageInx[], float value) {
        if (GFitsImageInx[0] == 1) {
            self->at(GFitsImageInx[1]) = value;
        }
        else if (GFitsImageInx[0] == 2) {
            self->at(GFitsImageInx[1], GFitsImageInx[2]) = value;
        }
        else if (GFitsImageInx[0] == 3) {
            self->at(GFitsImageInx[1], GFitsImageInx[2], GFitsImageInx[3]) =
                     value;
        }
        else if (GFitsImageInx[0] == 4) {
            self->at(GFitsImageInx[1], GFitsImageInx[2], GFitsImageInx[3],
                     GFitsImageInx[4]) = value;
        }
        else {
            throw GException::fits_wrong_image_operator("__setitem__(int)",
                                                        self->naxis(),
                                                        GFitsImageInx[0]);
        }
    }
    PyObject* pixels(void) {
        PyObject *list = PyList_New(self->npix());
        for (int i = 0; i < self->npix(); ++i) {
            PyList_SetItem(list, i, PyFloat_FromDouble((double)(static_cast<float*>(self->pixels())[i])));
        }
        return list;
    }
    GFitsImageFloat copy() {
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
