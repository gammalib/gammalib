/***************************************************************************
 *         GMWLException.hpp  - Multi-wavelength exception handler         *
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
 * @file GMWLException.hpp
 * @brief MWL exception handler interface definition.
 * @author J. Knodlseder
 */

#ifndef GMWLEXCEPTION_HPP
#define GMWLEXCEPTION_HPP

/* __ Includes ___________________________________________________________ */
#include <string>                             // string
#include <sstream>                            // ostringstream
#include <stdexcept>                          // exception


/***********************************************************************//**
 * @class GMWLExceptionHandler
 *
 * @brief Interface for multi-wavelength exception handler
 ***************************************************************************/
class GMWLExceptionHandler : public std::exception {
public:
    GMWLExceptionHandler(void) { }
    virtual ~GMWLExceptionHandler(void) throw() { }
    virtual const char* what() const throw();
protected:
    std::string m_origin;
    std::string m_message;
};


/***********************************************************************//**
 * @class GMWLException
 *
 * @brief Interface for multi-wavelength exceptions
 *
 * This object is thrown in case of a multi-wavelength exception.
 ***************************************************************************/
class GMWLException : public GMWLExceptionHandler {
public:
    // Exceptions
    class file_open_error : public GMWLExceptionHandler {
    public:
        file_open_error(std::string origin, std::string filename,
                        std::string message = "");
    };
    class bad_file_format : public GMWLExceptionHandler {
    public:
        bad_file_format(std::string origin, std::string message = "");
    };
};

#endif /* GMWLEXCEPTION_HPP */
