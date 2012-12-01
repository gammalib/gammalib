/***************************************************************************
 *                GCTAException.hpp  - CTA exception handler               *
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
 * @file GCTAException.hpp
 * @brief CTA exception handler interface definition.
 * @author Juergen Knoedlseder
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
    // General exceptions
    class file_open_error : public GCTAExceptionHandler {
    public:
        file_open_error(std::string origin, std::string filename,
                        std::string message = "");
    };

    // Event bin exceptions
    class no_member : public GCTAExceptionHandler {
    public:
        no_member(const std::string& origin, const std::string& message);
    };

    // Observation exceptions
    class no_response : public GCTAExceptionHandler {
    public:
        no_response(std::string origin, std::string message = "");
    };
    class no_pointing : public GCTAExceptionHandler {
    public:
        no_pointing(std::string origin, std::string message = "");
    };
    class no_sky : public GCTAExceptionHandler {
    public:
        no_sky(std::string origin, std::string message = "");
    };
    class no_ebds : public GCTAExceptionHandler {
    public:
        no_ebds(std::string origin, std::string message = "");
    };
    class no_gti : public GCTAExceptionHandler {
    public:
        no_gti(std::string origin, std::string message = "");
    };
    class no_energies : public GCTAExceptionHandler {
    public:
        no_energies(std::string origin, std::string message = "");
    };
    class no_dirs : public GCTAExceptionHandler {
    public:
        no_dirs(std::string origin, std::string message = "");
    };
    class bad_observation_type : public GCTAExceptionHandler {
    public:
        bad_observation_type(std::string origin, std::string message = "");
    };
    class bad_event_type : public GCTAExceptionHandler {
    public:
        bad_event_type(std::string origin, std::string message = "");
    };
    class bad_roi_type : public GCTAExceptionHandler {
    public:
        bad_roi_type(std::string origin, std::string message = "");
    };
    class bad_instdir_type : public GCTAExceptionHandler {
    public:
        bad_instdir_type(std::string origin, std::string message = "");
    };
    class bad_pointing_type : public GCTAExceptionHandler {
    public:
        bad_pointing_type(std::string origin, std::string message = "");
    };
    class bad_response_type : public GCTAExceptionHandler {
    public:
        bad_response_type(std::string origin, std::string message = "");
    };

    // Model exceptions
    class model_invalid_radial : public GCTAExceptionHandler {
    public:
        model_invalid_radial(std::string origin, std::string type,
                             std::string message = "");
    };

    // Response table exceptions
    class bad_rsp_table_dim : public GCTAExceptionHandler {
    public:
        bad_rsp_table_dim(std::string origin,
                          int         table_dim,
                          int         expected_dim,
                          std::string message = "");
    };
    class bad_rsp_table_format : public GCTAExceptionHandler {
    public:
        bad_rsp_table_format(std::string origin,
                             std::string message = "");
    };
};

#endif /* GCTAEXCEPTION_HPP */
