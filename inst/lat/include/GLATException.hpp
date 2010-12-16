/***************************************************************************
 *                GLATException.hpp  - LAT exception handler               *
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
 * @file GLATException.hpp
 * @brief LAT exception handler interface definition.
 * @author J. Knodlseder
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
    class no_response : public GLATExceptionHandler {
    public:
        no_response(std::string origin, std::string message = "");
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
};

#endif /* GLATEXCEPTION_HPP */
