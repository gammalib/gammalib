/***************************************************************************
 *                      GException.i - exception handler                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2012 by Juergen Knoedlseder                         *
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
 * @file GException.i
 * @brief Exception handler for Python interface
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here */
#include "GException.hpp"
%}
%include exception.i

%exception {
    try {
        $action
    }
    catch (const GException::out_of_range& e) {
        SWIG_exception(SWIG_IndexError, e.what());
    }
    catch (const GExceptionHandler& e) {
        SWIG_exception(SWIG_RuntimeError, e.what());
    }
    catch (const std::exception& e) {
        SWIG_exception(SWIG_RuntimeError, e.what());
    }
}
