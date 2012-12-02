/***************************************************************************
 *              GCOMException.hpp  - COMPTEL exception handler             *
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
 * @file GCOMException.hpp
 * @brief COMPTEL exception handler interface definition
 * @author J. Knoedlseder
 */

#ifndef GCOMEXCEPTION_HPP
#define GCOMEXCEPTION_HPP

/* __ Includes ___________________________________________________________ */
#include <string>                             // string
#include <sstream>                            // ostringstream
#include <stdexcept>                          // exception


/***********************************************************************//**
 * @class GCOMExceptionHandler
 *
 * @brief Interface for COMPTEL exception handler
 ***************************************************************************/
class GCOMExceptionHandler : public std::exception {
public:
    GCOMExceptionHandler(void) { }
    virtual ~GCOMExceptionHandler(void) throw() { }
    virtual const char* what() const throw();
protected:
    std::string m_origin;
    std::string m_message;
};


/***********************************************************************//**
 * @class GCOMException
 *
 * @brief Interface for COMPTEL exceptions
 *
 * This object is thrown in case of a COMPTEL exception.
 ***************************************************************************/
class GCOMException : public GCOMExceptionHandler {
public:
    // General exceptions
    class file_open_error : public GCOMExceptionHandler {
    public:
        file_open_error(const std::string& origin,
                        const std::string& filename,
                        const std::string& message = "");
    };

    // Event cube exceptions
    class no_sky : public GCOMExceptionHandler {
    public:
        no_sky(const std::string& origin,
               const std::string& message = "");
    };
    class no_ebds : public GCOMExceptionHandler {
    public:
        no_ebds(const std::string& origin,
                const std::string& message = "");
    };
    class no_gti : public GCOMExceptionHandler {
    public:
        no_gti(const std::string& origin,
               const std::string& message = "");
    };
    class no_dirs : public GCOMExceptionHandler {
    public:
        no_dirs(const std::string& origin,
                const std::string& message = "");
    };

    // Response exceptions
    class bad_observation_type : public GCOMExceptionHandler {
    public:
        bad_observation_type(const std::string& origin,
                             const std::string& message = "");
    };
    class bad_event_type : public GCOMExceptionHandler {
    public:
        bad_event_type(const std::string& origin,
                       const std::string& message = "");
    };
};

#endif /* GCOMEXCEPTION_HPP */
