/***************************************************************************
 *                 typemap_slices.i - Typemaps for slicing                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017 by Juergen Knoedlseder                              *
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
 * @file typemap_slices.i
 * @brief Typemaps for slicing
 * @author Juergen Knoedlseder
 */

/***********************************************************************//**
 * @brief Get slice information from Python object
 *
 * @param[in] slice Python object.
 * @param[in] length Sequence length.
 * @param[out] start Start index.
 * @param[out] stop Stop index.
 * @param[out] step Steo size.
 * @return Return code (0: success, -1: failure)
 ***************************************************************************/
%inline %{
int PythonSlice_GetIndices(PyObject *  slice,
                           Py_ssize_t  length,
                           Py_ssize_t* start,
                           Py_ssize_t* stop,
                           Py_ssize_t* step) {
    #if PY_VERSION_HEX >= 0x03020000
    int rc = PySlice_GetIndices(slice, length, start, stop, step);
    #else
    int rc = PySlice_GetIndices((PySliceObject*)slice, length, start, stop, step);
    #endif
    return rc;
}
%}
