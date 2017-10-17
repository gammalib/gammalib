/***************************************************************************
 *                         GSkyMap.i - Sky map class                       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2017 by Juergen Knoedlseder                         *
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
 * @file GSkyMap.i
 * @brief Sky map class SWIG file.
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GSkyMap.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @brief Tuple to index conversion to provide pixel access.
 *
 * The following function provides conversion between a Python tuple and
 * an integer array. This allows skymap pixel access via tuples, such as in
 * a[3,5,10] = 10.0 or c = a[2,9].
 ***************************************************************************/
%{
static int skymap_tuple(PyObject *input, int *ptr) {
    if (PySequence_Check(input)) {
        int size = PyObject_Length(input);
        if (size > 2) {
            PyErr_SetString(PyExc_ValueError,"Too many arguments in tuple");
            return 0;
        }
        ptr[0] = size;
        for (int i = 0; i < size; i++) {
            PyObject *o = PySequence_GetItem(input,i);
            if (!PyInt_Check(o)) {
                Py_XDECREF(o);
                PyErr_SetString(PyExc_ValueError,"Expecting a tuple of integers");
                return 0;
            }
            ptr[i+1] = (int)PyInt_AsLong(o);
            Py_DECREF(o);
        }
        return 1;
    }
    else {
        ptr[0] = 1;
        if (!PyInt_Check(input)) {
            PyErr_SetString(PyExc_ValueError,"Expecting an integer");
            return 0;
        }
        ptr[1] = (int)PyInt_AsLong(input);
        return 1;       
    }
}
%}

// This is the typemap that makes use of the function defined above
%typemap(in) int GSkyMapInx[ANY](int temp[3]) {
   if (!skymap_tuple($input,temp)) {
      return NULL;
   }
   $1 = &temp[0];
}

// This typecheck verifies that all arguments are integers. The typecheck
// is needed for using "int GSkyMapInx" in overloaded methods.
%typemap(typecheck) int GSkyMapInx[ANY] {
    $1 = 1;
    if (PySequence_Check($input)) {
        int size = PyObject_Length($input);
        for (int i = 0; i < size; i++) {
            PyObject *o = PySequence_GetItem($input,i);
            if (!PyInt_Check(o)) {
                $1 = 0;
                break;
            }
        }
    }
    else {
        if (!PyInt_Check($input)) {
            $1 = 0;
        }
    }
}


/***********************************************************************//**
 * @class GSkyMap
 *
 * @brief GSkyMap class interface definition
 ***************************************************************************/
class GSkyMap : public GBase {

public:
    // Constructors and destructors
    GSkyMap(void);
    explicit GSkyMap(const GFilename& filename);
    explicit GSkyMap(const GFitsHDU& hdu);
    GSkyMap(const std::string& coords,
            const int&         nside,
            const std::string& order,
            const int&         nmaps = 1);
    GSkyMap(const std::string& wcs,
            const std::string& coords,
            const double&      x,
            const double&      y,
            const double&      dx,
            const double&      dy,
            const int&         nx,
            const int&         ny,
            const int&         nmaps = 1);
    GSkyMap(const GSkyMap& map);
    virtual ~GSkyMap(void);

    // Operators
    GSkyMap&      operator+=(const GSkyMap& map);
    GSkyMap&      operator+=(const double& value);
    GSkyMap&      operator-=(const GSkyMap& map);
    GSkyMap&      operator-=(const double& value);
    GSkyMap&      operator*=(const GSkyMap& map);
    GSkyMap&      operator*=(const double& factor);
    GSkyMap       operator+(const GSkyMap& map) const;
    GSkyMap       operator-(const GSkyMap& map) const;
    GSkyMap       operator*(const GSkyMap& map) const;
    GSkyMap       operator/(const GSkyMap& map) const;
    double        operator()(const int& index, const int& map = 0);
    double        operator()(const GSkyPixel& pixel, const int& map = 0);
    double        operator()(const GSkyDir& dir, const int& map = 0) const;

    // Methods
    void                    clear(void);
    GSkyMap*                clone(void) const;
    std::string             classname(void) const;
    bool                    is_empty(void) const;
    const int&              npix(void) const;
    const int&              nx(void) const;
    const int&              ny(void) const;
    const int&              nmaps(void) const;
    void                    nmaps(const int& nmaps);
    const std::vector<int>& shape(void) const;
    void                    shape(const int& s1);
    void                    shape(const int& s1, const int& s2);
    void                    shape(const int& s1, const int& s2, const int& s3);
    void                    shape(const std::vector<int>& shape);
    int                     ndim(void) const;
    GSkyPixel               inx2pix(const int& index) const;
    GSkyDir                 inx2dir(const int& index) const;
    GSkyDir                 pix2dir(const GSkyPixel& pixel) const;
    int                     pix2inx(const GSkyPixel& pixel) const;
    int                     dir2inx(const GSkyDir& dir) const;
    GSkyPixel               dir2pix(const GSkyDir& dir) const;
    double                  flux(const int& index, const int& map = 0) const;
    double                  flux(const GSkyPixel& pixel, const int& map = 0) const;
    double                  solidangle(const int& index) const;
    double                  solidangle(const GSkyPixel& pixel) const;
    bool                    contains(const GSkyDir& dir) const;
    bool                    contains(const GSkyPixel& pixel) const;
    bool                    overlaps(const GSkyRegion& region) const;
    void                    smooth(const std::string& kernel, const double& par);
    const GSkyProjection*   projection(void) const;
    void                    projection(const GSkyProjection& proj);
    const double*           pixels(void) const;
    GSkyMap                 extract(const int& map, const int& nmaps = 1) const;
    void                    stack_maps(void);
    void                    load(const GFilename& filename);
    void                    save(const GFilename& filename,
                                 const bool&      clobber = false) const;
    void                    read(const GFitsHDU& hdu);
    GFitsHDU*               write(GFits& file,
                                  const std::string& extname = "") const;
    void                    publish(const std::string& name = "") const;
};


/***********************************************************************//**
 * @brief GSkyMap class extension
 *
 * @todo Implement __getitem__ and __setitem__ methods for GSkyPixel and
 * GSkyDir
 ***************************************************************************/
%extend GSkyMap {
    double __getitem__(int GSkyMapInx[]) {
        if (GSkyMapInx[1] < 0 || GSkyMapInx[1] >= self->npix()) {
            throw GException::out_of_range("__getitem__(int)", "Sky map index",
                                           GSkyMapInx[1], self->npix());
        }
        if (GSkyMapInx[0] == 1) {
            return (*self)(GSkyMapInx[1]);
        }
        else {
            if (GSkyMapInx[2] >= 0 && GSkyMapInx[2] < self->nmaps()) {
                return (*self)(GSkyMapInx[1], GSkyMapInx[2]);
            }
            else {
                throw GException::out_of_range("__getitem__(int)", "Sky map number",
                                               GSkyMapInx[2], self->nmaps());
            }
        }
    }
    /*
    double __getitem__(const GSkyPixel& pixel) {
        return (*self)(pixel);
    }
    */
    void __setitem__(int GSkyMapInx[], double value) {
        if (GSkyMapInx[1] < 0 || GSkyMapInx[1] >= self->npix()) {
            throw GException::out_of_range("__setitem__(int)", "Sky map index",
                                           GSkyMapInx[1], self->npix());
        }
        if (GSkyMapInx[0] == 1) {
            (*self)(GSkyMapInx[1]) = value;
        }
        else {
            if (GSkyMapInx[2] >= 0 && GSkyMapInx[2] < self->nmaps()) {
                (*self)(GSkyMapInx[1], GSkyMapInx[2]) = value;
            }
            else {
                throw GException::out_of_range("__setitem__(int)", "Sky map number",
                                               GSkyMapInx[2], self->nmaps());
            }
        }
    }
    /*
    void __setitem__(const GSkyPixel& pixel, double value) {
        (*self)(pixel) = value;
    }
    */
    GSkyMap copy() {
        return (*self);
    }
    GSkyMap sqrt() {
        return sqrt(*self);
    }
    GSkyMap log() {
        return log(*self);
    }
    GSkyMap log10() {
        return log10(*self);
    }
    GSkyMap abs() {
        return abs(*self);
    }
    GSkyMap sign() {
        return sign(*self);
    }
    // Python 2.x operator/=
    GSkyMap __idiv__(const GSkyMap& map) {
        self->operator/=(map);
        return (*self);
    }
    // Python 3.x operator/=
    GSkyMap __itruediv__(const GSkyMap& map) {
        self->operator/=(map);
        return (*self);
    }
    // Python 2.x operator/=
    GSkyMap __idiv__(const double& factor) {
        self->operator/=(factor);
        return (*self);
    }
    // Python 3.x operator/=
    GSkyMap __itruediv__(const double& factor) {
        self->operator/=(factor);
        return (*self);
    }
};
