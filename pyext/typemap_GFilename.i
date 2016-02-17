/***************************************************************************
 *                  typemap_GFilename.i - Filename typemap                 *
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
 * @file typemap_GFilename.i
 * @brief GFilename SWIG typemap
 * @author Juergen Knoedlseder
 */

/***********************************************************************//**
 * @brief Input typemap for GFilename
 *
 * This typemap allows to pass a string for a GFilename& argument. In that
 * way the user can simply pass a string to any function or method that has
 * a GFilename argument, and the typemap transparently converts a string
 * into a GFilename object.
 *
 * Note that the $1_descriptor special variable will check for the GFilename
 * pointer.
 ***************************************************************************/
%typemap(in) GFilename& (GFilename temp) {
    if (PyString_Check($input)) {
        temp = GFilename(std::string(PyString_AsString($input)));
        $1 = &temp;
    }
    else if (PyUnicode_Check($input)) {
        temp = GFilename(std::string(PyBytes_AS_STRING(PyUnicode_AsEncodedString($input, "utf-8", "Error ~"))));
        $1 = &temp;
    }
    else {
        void *filename_argp1 = 0;
        if (SWIG_IsOK(SWIG_ConvertPtr($input, &filename_argp1, $1_descriptor, 0))) {
            $1 = reinterpret_cast<GFilename*>(filename_argp1);
        }
        else {
            SWIG_exception(SWIG_TypeError, "GFilename expected");
        }
    }
}


/***********************************************************************//**
 * @brief Typecheck typemap for GFilename
 *
 * This typemap allows to overload methods or constructors and to pass a
 * string for a GFilename& argument. We take here as precedence a typecheck
 * for a double precision value. A precedence is needed as otherwise a
 * warning will be issued.
 ***************************************************************************/
%typecheck(SWIG_TYPECHECK_DOUBLE) GFilename& {
    if (PyString_Check($input) ||
        SWIG_CheckState(SWIG_ConvertPtr($input, 0, $1_descriptor, 0))) {
        $1 = 1;
    }
    else if (PyUnicode_Check($input) ||
        SWIG_CheckState(SWIG_ConvertPtr($input, 0, $1_descriptor, 0))) {
        $1 = 1;
    }
    else {
        $1 = 0;
    }
}
