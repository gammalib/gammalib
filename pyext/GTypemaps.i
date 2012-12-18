/***************************************************************************
 *          GTypemaps.i  -  Typemaps for GammaLib Python interface         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012 by Juergen Knoedlseder                              *
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
 * @file GTypemaps.i
 * @brief Provides typemaps for GammaLib
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GFitsHDU.hpp"
#include "GFitsImage.hpp"
#include "GFitsImageByte.hpp"
#include "GFitsImageSByte.hpp"
#include "GFitsImageUShort.hpp"
#include "GFitsImageShort.hpp"
#include "GFitsImageULong.hpp"
#include "GFitsImageLong.hpp"
#include "GFitsImageLongLong.hpp"
#include "GFitsImageFloat.hpp"
#include "GFitsImageDouble.hpp"
#include "GFitsTable.hpp"
#include "GFitsAsciiTable.hpp"
#include "GFitsBinTable.hpp"
#include "GModelSky.hpp"
#include "GModelDiffuseSource.hpp"
#include "GModelExtendedSource.hpp"
#include "GModelPointSource.hpp"
#include "GModelData.hpp"
#include "GModelSpatial.hpp"
#include "GModelRadial.hpp"
#include "GModelRadialDisk.hpp"
#include "GModelRadialGauss.hpp"
#include "GModelRadialShell.hpp"
#include "GModelSpatialConst.hpp"
#include "GModelSpatialCube.hpp"
#include "GModelSpatialMap.hpp"
#include "GModelSpatialPtsrc.hpp"
#include "GModelSpectral.hpp"
#include "GModelSpectralConst.hpp"
#include "GModelSpectralFunc.hpp"
#include "GModelSpectralNodes.hpp"
#include "GModelSpectralPlaw.hpp"
#include "GModelSpectralPlaw2.hpp"
#include "GModelTemporal.hpp"
#include "GModelTemporalConst.hpp"
%}


/***********************************************************************//**
 * @brief 2D tuple to index conversion
 *
 * The following typemap provides conversion between a 2-dimensional Python
 * tuple and an integer array. This allows index access via tuples, such as
 * in a[3,5] = 10.0 or c = a[2,9]. A typecheck typemap is provided to allow
 * overloading.
 ***************************************************************************/
%typemap(in) int GTuple[2] (int temp[2]) {
    if (!PySequence_Check($input)) {
        PyErr_SetString(PyExc_ValueError,"Expected a sequence");
        return NULL;
    }
    if (PySequence_Length($input) != 2) {
        PyErr_SetString(PyExc_ValueError,"Size mismatch. Expected 2 elements");
        return NULL;
    }
    for (int i = 0; i < 2; ++i) {
        PyObject *o = PySequence_GetItem($input,i);
        if (PyInt_Check(o)) {
            temp[i] = (int)PyInt_AsLong(o);
        } 
        else {
            PyErr_SetString(PyExc_ValueError,"Indices must be integers");      
            return NULL;
        }
    }
    $1 = temp;
}
%typemap(typecheck) int GTuple[2] {
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
 * @brief Tuple to index conversion using variable dimensions (1,...,4)
 *
 * The following function provides conversion between a Python tuple and
 * an integer array. This allows index access via tuples, such as in
 * a[3,5,10] = 10.0 or c = a[2,9].
 ***************************************************************************/
%{
static int var_tuple_to_index(PyObject *input, int *ptr) {
    if (PySequence_Check(input)) {
        int size = PyObject_Length(input);
        if (size > 4) {
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
%typemap(in) int GTuple[ANY] (int temp[5]) {
   if (!var_tuple_to_index($input,temp)) {
      return NULL;
   }
   $1 = &temp[0];
}
%typemap(typecheck) int GTuple[ANY] {
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
 * @brief GFitsHDU* output conversion
 *
 * This typemap implements an automatic cast of a GFitsHDU pointer to
 * the relevant derived class pointer.
 ***************************************************************************/
%typemap(out) GFitsHDU* {
    if (dynamic_cast<GFitsImage*>($1) != NULL) {
        if (dynamic_cast<GFitsImageByte*>($1) != NULL) {
            $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsImageByte, 0 |  0 );
        }
        else if (dynamic_cast<GFitsImageDouble*>($1) != NULL) {
            $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsImageDouble, 0 |  0 );
        }
        else if (dynamic_cast<GFitsImageFloat*>($1) != NULL) {
            $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsImageFloat, 0 |  0 );
        }
        else if (dynamic_cast<GFitsImageLong*>($1) != NULL) {
            $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsImageFloat, 0 |  0 );
        }
        else if (dynamic_cast<GFitsImageLongLong*>($1) != NULL) {
            $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsImageLongLong, 0 |  0 );
        }
        else if (dynamic_cast<GFitsImageSByte*>($1) != NULL) {
            $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsImageSByte, 0 |  0 );
        }
        else if (dynamic_cast<GFitsImageShort*>($1) != NULL) {
            $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsImageShort, 0 |  0 );
        }
        else if (dynamic_cast<GFitsImageULong*>($1) != NULL) {
            $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsImageULong, 0 |  0 );
        }
        else if (dynamic_cast<GFitsImageUShort*>($1) != NULL) {
            $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsImageUShort, 0 |  0 );
        }
        else {
            $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsImage, 0 |  0 );
        }
    }
    else if (dynamic_cast<GFitsTable*>($1) != NULL) {
        if (dynamic_cast<GFitsBinTable*>($1) != NULL) {
            $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsBinTable, 0 |  0 );
        }
        else if (dynamic_cast<GFitsAsciiTable*>($1) != NULL) {
            $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsAsciiTable, 0 |  0 );
        }
        else {
            $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsTable, 0 |  0 );
        }
    }
    else {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsHDU, 0 |  0 );
    }
}


/***********************************************************************//**
 * @brief GFitsImage* output conversion
 *
 * This typemap implements an automatic cast of a GFitsImage pointer to
 * the relevant derived class pointer.
 ***************************************************************************/
%typemap(out) GFitsImage* {
    if (dynamic_cast<GFitsImageByte*>($1) != NULL) {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsImageByte, 0 |  0 );
    }
    else if (dynamic_cast<GFitsImageDouble*>($1) != NULL) {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsImageDouble, 0 |  0 );
    }
    else if (dynamic_cast<GFitsImageFloat*>($1) != NULL) {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsImageFloat, 0 |  0 );
    }
    else if (dynamic_cast<GFitsImageLong*>($1) != NULL) {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsImageFloat, 0 |  0 );
    }
    else if (dynamic_cast<GFitsImageLongLong*>($1) != NULL) {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsImageLongLong, 0 |  0 );
    }
    else if (dynamic_cast<GFitsImageSByte*>($1) != NULL) {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsImageSByte, 0 |  0 );
    }
    else if (dynamic_cast<GFitsImageShort*>($1) != NULL) {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsImageShort, 0 |  0 );
    }
    else if (dynamic_cast<GFitsImageULong*>($1) != NULL) {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsImageULong, 0 |  0 );
    }
    else if (dynamic_cast<GFitsImageUShort*>($1) != NULL) {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsImageUShort, 0 |  0 );
    }
    else {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsImage, 0 |  0 );
    }
}


/***********************************************************************//**
 * @brief GFitsTable* output conversion
 *
 * This typemap implements an automatic cast of a GFitsTable pointer to
 * the relevant derived class pointer.
 ***************************************************************************/
%typemap(out) GFitsTable* {
    if (dynamic_cast<GFitsBinTable*>($1) != NULL) {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsBinTable, 0 |  0 );
    }
    else if (dynamic_cast<GFitsAsciiTable*>($1) != NULL) {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsAsciiTable, 0 |  0 );
    }
    else {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GFitsTable, 0 |  0 );
    }
}


/***********************************************************************//**
 * @brief GModel& output conversion
 *
 * This typemap implements an automatic cast of a GModel reference to
 * the relevant derived class reference.
 ***************************************************************************/
%typemap(out) GModel* {
    if (dynamic_cast<GModelSky*>($1) != NULL) {
        if (dynamic_cast<GModelDiffuseSource*>($1) != NULL) {
            $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GModelDiffuseSource, 0 |  0 );
        }
        else if (dynamic_cast<GModelExtendedSource*>($1) != NULL) {
            $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GModelExtendedSource, 0 |  0 );
        }
        else if (dynamic_cast<GModelPointSource*>($1) != NULL) {
            $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GModelPointSource, 0 |  0 );
        }
        else {
            $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GModelSky, 0 |  0 );
        }
    }
    else if (dynamic_cast<GModelData*>($1) != NULL) {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GModelData, 0 |  0 );
    }
    else {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GModel, 0 |  0 );
    }
}


/***********************************************************************//**
 * @brief GModelSpatial* output conversion
 *
 * This typemap implements an automatic cast of a GModelSpatial pointer to
 * the relevant derived class pointer.
 ***************************************************************************/
%typemap(out) GModelSpatial* {
    if (dynamic_cast<GModelRadial*>($1) != NULL) {
        if (dynamic_cast<GModelRadialDisk*>($1) != NULL) {
            $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GModelRadialDisk, 0 |  0 );
        }
        else if (dynamic_cast<GModelRadialGauss*>($1) != NULL) {
            $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GModelRadialGauss, 0 |  0 );
        }
        else if (dynamic_cast<GModelRadialShell*>($1) != NULL) {
            $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GModelRadialShell, 0 |  0 );
        }
        else {
            $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GModelRadial, 0 |  0 );
        }
    }
    else if (dynamic_cast<GModelSpatialConst*>($1) != NULL) {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GModelSpatialConst, 0 |  0 );
    }
    else if (dynamic_cast<GModelSpatialCube*>($1) != NULL) {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GModelSpatialCube, 0 |  0 );
    }
    else if (dynamic_cast<GModelSpatialMap*>($1) != NULL) {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GModelSpatialMap, 0 |  0 );
    }
    else if (dynamic_cast<GModelSpatialPtsrc*>($1) != NULL) {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GModelSpatialPtsrc, 0 |  0 );
    }
    else {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GModelSpatial, 0 |  0 );
    }
}


/***********************************************************************//**
 * @brief GModelSpectral* output conversion
 *
 * This typemap implements an automatic cast of a GModelSpectral pointer to
 * the relevant derived class pointer.
 ***************************************************************************/
%typemap(out) GModelSpectral* {
    if (dynamic_cast<GModelSpectralConst*>($1) != NULL) {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GModelSpectralConst, 0 |  0 );
    }
    else if (dynamic_cast<GModelSpectralFunc*>($1) != NULL) {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GModelSpectralFunc, 0 |  0 );
    }
    else if (dynamic_cast<GModelSpectralNodes*>($1) != NULL) {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GModelSpectralNodes, 0 |  0 );
    }
    else if (dynamic_cast<GModelSpectralPlaw*>($1) != NULL) {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GModelSpectralPlaw, 0 |  0 );
    }
    else if (dynamic_cast<GModelSpectralPlaw2*>($1) != NULL) {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GModelSpectralPlaw2, 0 |  0 );
    }
    else if (dynamic_cast<GModelSpectralLogParabola*>($1) != NULL) {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GModelSpectralLogParabola, 0 |  0 );
    }
    else {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GModelSpectral, 0 |  0 );
    }
}


/***********************************************************************//**
 * @brief GModelTemporal* output conversion
 *
 * This typemap implements an automatic cast of a GModelTemporal pointer to
 * the relevant derived class pointer.
 ***************************************************************************/
%typemap(out) GModelTemporal* {
    if (dynamic_cast<GModelTemporalConst*>($1) != NULL) {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GModelTemporalConst, 0 |  0 );
    }
    else {
        $result = SWIG_NewPointerObj(SWIG_as_voidptr($1), SWIGTYPE_p_GModelTemporal, 0 |  0 );
    }
}
