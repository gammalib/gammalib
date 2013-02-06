/***************************************************************************
 *                GXXXException.hpp  - XXX exception handler               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013 by Juergen Knoedlseder                              *
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
 * @file GXXXException.hpp
 * @brief XXX exception handler interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GXXXEXCEPTION_HPP
#define GXXXEXCEPTION_HPP

/* __ Includes ___________________________________________________________ */
#include <string>                             // string
#include <sstream>                            // ostringstream
#include <stdexcept>                          // exception


/***********************************************************************//**
 * @class GXXXExceptionHandler
 *
 * @brief Interface for XXX exception handler
 ***************************************************************************/
class GXXXExceptionHandler : public std::exception {
public:
    GXXXExceptionHandler(void) { }
    virtual ~GXXXExceptionHandler(void) throw() { }
    virtual const char* what() const throw();
protected:
    std::string m_origin;
    std::string m_message;
};


/***********************************************************************//**
 * @class GXXXException
 *
 * @brief Interface for XXX exceptions
 *
 * This object is thrown in case of a XXX exception.
 ***************************************************************************/
class GXXXException : public GXXXExceptionHandler {
public:
};

#endif /* GXXXEXCEPTION_HPP */
