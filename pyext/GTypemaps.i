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
