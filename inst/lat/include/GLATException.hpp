/***************************************************************************
 *                GLATException.hpp  - LAT exception handler               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2012 by Juergen Knoedlseder                         *
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
 * @file GLATException.hpp
 * @brief LAT exception handler interface definition.
 * @author Juergen Knoedlseder
 */

#ifndef GLATEXCEPTION_HPP
#define GLATEXCEPTION_HPP

/* __ Includes ___________________________________________________________ */
#include <string>                             // string
#include <sstream>                            // ostringstream
#include <stdexcept>                          // exception


/***********************************************************************//**
 * @class GLATExceptionHandler
 *
 * @brief Interface for LAT exception handler.
 ***************************************************************************/
class GLATExceptionHandler : public std::exception {
public:
    GLATExceptionHandler(void) { }
    virtual ~GLATExceptionHandler(void) throw() { }
    virtual const char* what() const throw();
protected:
    std::string m_origin;
    std::string m_message;
};


/***********************************************************************//**
 * @class GLATException
 *
 * @brief Interface for LAT exceptions.
 *
 * This is the object that is thrown in case of a LAT exception.
 ***************************************************************************/
class GLATException : public GLATExceptionHandler {
public:
    // Exceptions
    class file_open_error : public GLATExceptionHandler {
    public:
        file_open_error(std::string origin, std::string filename,
                        std::string message = "");
    };
    class no_member : public GLATExceptionHandler {
    public:
        no_member(const std::string& origin, const std::string& message);
    };
    class no_sky : public GLATExceptionHandler {
    public:
        no_sky(std::string origin, std::string message = "");
    };
    class no_ebds : public GLATExceptionHandler {
    public:
        no_ebds(std::string origin, std::string message = "");
    };
    class no_gti : public GLATExceptionHandler {
    public:
        no_gti(std::string origin, std::string message = "");
    };
    class no_ltcube : public GLATExceptionHandler {
    public:
        no_ltcube(std::string origin, std::string message = "");
    };
    class no_energies : public GLATExceptionHandler {
    public:
        no_energies(std::string origin, std::string message = "");
    };
    class no_dirs : public GLATExceptionHandler {
    public:
        no_dirs(std::string origin, std::string message = "");
    };
    class bad_roi_type : public GLATExceptionHandler {
    public:
        bad_roi_type(std::string origin, std::string message = "");
    };
    class bad_instdir_type : public GLATExceptionHandler {
    public:
        bad_instdir_type(std::string origin, std::string message = "");
    };
    class invalid_response : public GLATExceptionHandler {
    public:
        invalid_response(std::string origin, std::string message = "");
    };
    class wcs_incompatible : public GLATExceptionHandler {
    public:
        wcs_incompatible(std::string origin, std::string name,
                         std::string message = "");
    };
    class diffuse_not_found : public GLATExceptionHandler {
    public:
        diffuse_not_found(std::string origin, std::string name,
                          std::string message = "");
    };
    class inconsistent_response : public GLATExceptionHandler {
    public:
        inconsistent_response(std::string origin, int size, int expect,
                              std::string message = "");
    };
    class bad_response_type : public GLATExceptionHandler {
    public:
        bad_response_type(std::string origin, std::string message = "");
    };
};

#endif /* GLATEXCEPTION_HPP */
