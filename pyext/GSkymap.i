/***************************************************************************
 *                GSkymap.i  -  Sky map class SWIG definition              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GSkymap.i
 * @brief GSkymap class SWIG file.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GSkymap.hpp"
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
%typemap(in) int GSkymapInx[ANY](int temp[3]) {
   if (!skymap_tuple($input,temp)) {
      return NULL;
   }
   $1 = &temp[0];
}


/***********************************************************************//**
 * @class GSkymap
 *
 * @brief GSkymap class interface defintion
 ***************************************************************************/
class GSkymap {
public:
    // Constructors and destructors
    GSkymap(void);
    explicit GSkymap(const std::string& filename);
    explicit GSkymap(const std::string& wcs, const std::string& coords,
                     const int& nside, const std::string& order,
                     const int nmaps = 1);
    explicit GSkymap(const std::string& wcs, const std::string& coords,
                     double const& x, double const& y,
                     double const& dx, double const& dy,
                     const int& nx, const int& ny, const int nmaps = 1);
    GSkymap(const GSkymap& map);
    virtual ~GSkymap(void);

    // Methods
    GSkyDir   pix2dir(const int& pix);
    int       dir2pix(GSkyDir dir) const;
    double    omega(const int& pix) const;
    GSkyDir   xy2dir(const GSkyPixel& pix);
    GSkyPixel dir2xy(GSkyDir dir) const;
    void      load(const std::string& filename);
    void      save(const std::string& filename, bool clobber = false) const;
    void      read(const GFitsHDU* hdu);
    void      write(GFits* file) const;
    int       npix(void) const;
    int       nx(void) const;
    int       ny(void) const;
    int       nmaps(void) const;
    GWcs*     wcs(void) const { return m_wcs; }
    double*   pixels(void) const { return m_pixels; }
};


/***********************************************************************//**
 * @brief GSkymap class extension
 ***************************************************************************/
%extend GSkymap {
    char *__str__() {
        static std::string result = self->print();
        return ((char*)result.c_str());
    }
    double __getitem__(int GSkymapInx[]) {
        if (GSkymapInx[0] == 1)
            return (*self)(GSkymapInx[1]);
        else
            return (*self)(GSkymapInx[1], GSkymapInx[2]);
    }
    /*
    double __getitem__(const GSkyPixel& pixel) {
        return (*self)(pixel);
    }
    */
    void __setitem__(int GSkymapInx[], double value) {
        if (GSkymapInx[0] == 1)
            (*self)(GSkymapInx[1]) = value;
        else
            (*self)(GSkymapInx[1], GSkymapInx[2]) = value;
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
