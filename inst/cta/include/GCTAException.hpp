/***************************************************************************
 *                GCTAException.hpp  - CTA exception handler               *
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
 * @file GCTAException.hpp
 * @brief CTA exception handler interface definition.
 * @author J. Knodlseder
 */

#ifndef GCTAEXCEPTION_HPP
#define GCTAEXCEPTION_HPP

/* __ Includes ___________________________________________________________ */
#include <string>                             // string
#include <sstream>                            // ostringstream
#include <stdexcept>                          // exception


/***********************************************************************//**
 * @class GCTAExceptionHandler
 *
 * @brief Interface for CTA exception handler.
 ***************************************************************************/
class GCTAExceptionHandler : public std::exception {
public:
    GCTAExceptionHandler(void) { }
    virtual ~GCTAExceptionHandler(void) throw() { }
    virtual const char* what() const throw();
protected:
    std::string m_origin;
    std::string m_message;
};


/***********************************************************************//**
 * @class GCTAException
 *
 * @brief Interface for CTA exceptions.
 *
 * This is the object that it thrown in case of a CTA exception.
 ***************************************************************************/
class GCTAException : public GCTAExceptionHandler {
public:
    // Exceptions
    class file_open_error : public GCTAExceptionHandler {
    public:
        file_open_error(std::string origin, std::string filename,
                        std::string message = "");
    };
    class response_not_set : public GCTAExceptionHandler {
    public:
        response_not_set(std::string origin,  std::string message = "");
    };
};

#endif /* GCTAEXCEPTION_HPP */
