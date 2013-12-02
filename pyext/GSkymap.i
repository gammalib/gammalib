/***************************************************************************
 *                         GSkymap.i - Sky map class                       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2013 by Juergen Knoedlseder                         *
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
 * @file GSkymap.i
 * @brief Sky map class SWIG file.
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GSkymap.hpp"
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
%typemap(in) int GSkymapInx[ANY](int temp[3]) {
   if (!skymap_tuple($input,temp)) {
      return NULL;
   }
   $1 = &temp[0];
}

// This typecheck verifies that all arguments are integers. The typecheck
// is needed for using "int GSkymapInx" in overloaded methods.
%typemap(typecheck) int GSkymapInx[ANY] {
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
 * @class GSkymap
 *
 * @brief GSkymap class interface defintion
 ***************************************************************************/
class GSkymap : public GBase {

public:
    // Constructors and destructors
    GSkymap(void);
    explicit GSkymap(const std::string& filename);
    explicit GSkymap(const std::string& proj, const std::string& coords,
                     const int& nside, const std::string& order,
                     const int nmaps = 1);
    explicit GSkymap(const std::string& proj, const std::string& coords,
                     double const& x, double const& y,
                     double const& dx, double const& dy,
                     const int& nx, const int& ny, const int nmaps = 1);
    GSkymap(const GSkymap& map);
    virtual ~GSkymap(void);

    // 1D pixel methods
    GSkyDir       pix2dir(const int& pix) const;
    int           dir2pix(const GSkyDir& dir) const;
    double        omega(const int& pix) const;

    // 2D pixel methods
    GSkyDir       xy2dir(const GSkyPixel& pix) const;
    GSkyPixel     dir2xy(const GSkyDir& dir) const;
    double        omega(const GSkyPixel& pix) const;

    // Methods
    void            clear(void);
    GSkymap*        clone(void) const;
    void            load(const std::string& filename);
    void            save(const std::string& filename, bool clobber = false) const;
    void            read(const GFitsHDU* hdu);
    void            write(GFits* file) const;
    int             npix(void) const;
    int             nx(void) const;
    int             ny(void) const;
    int             nmaps(void) const;
    int             xy2pix(const GSkyPixel& pix) const;
    GSkyPixel       pix2xy(const int& pix) const;
    GSkyProjection* projection(void) const;
    void            projection(const GSkyProjection& proj);
    double*         pixels(void) const;
    bool            contains(const GSkyDir& dir) const;
    bool            contains(const GSkyPixel& pixel) const;
};


/***********************************************************************//**
 * @brief GSkymap class extension
 *
 * @todo Implement __getitem__ and __setitem__ methods for GSkyPixel and
 * GSkyDir
 ***************************************************************************/
%extend GSkymap {
    double __getitem__(int GSkymapInx[]) {
        if (GSkymapInx[0] == 1) {
            return (*self)(GSkymapInx[1]);
        }
        else {
            return (*self)(GSkymapInx[1], GSkymapInx[2]);
        }
    }
    /*
    double __getitem__(const GSkyPixel& pixel) {
        return (*self)(pixel);
    }
    */
    void __setitem__(int GSkymapInx[], double value) {
        if (GSkymapInx[0] == 1) {
            (*self)(GSkymapInx[1]) = value;
        }
        else {
            (*self)(GSkymapInx[1], GSkymapInx[2]) = value;
        }
    }
    /*
    void __setitem__(const GSkyPixel& pixel, double value) {
        (*self)(pixel) = value;
    }
    */
    GSkymap copy() {
        return (*self);
    }
};
