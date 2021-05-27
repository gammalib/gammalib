/***************************************************************************
 *                    GException.hpp - Exception handler                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2006-2021 by Juergen Knoedlseder                         *
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
 * @file GException.hpp
 * @brief Exception handler interface definition.
 * @author Juergen Knoedlseder
 */

#ifndef GEXCEPTION_HPP
#define GEXCEPTION_HPP

/* __ Includes ___________________________________________________________ */
#include <string>                             // string
#include <stdexcept>                          // exception

/* __ Forward declarations ______________________________________________ */
class GEnergy;

/* __ Prototypes ________________________________________________________ */
namespace gammalib {
    void check_energy_interval(const std::string& origin,
                               const GEnergy&     emin,
                               const GEnergy&     emax);
    void check_prj_x2s_status(const std::string& origin,
                              const int&         status,
                              const int&         number);
    void check_prj_s2x_status(const std::string& origin,
                              const int&         status,
                              const int&         number);
}


/***********************************************************************//**
 * @class GExceptionHandler
 *
 * @brief Interface for exception handler.
 ***************************************************************************/
class GExceptionHandler : public std::exception {
public:
    GExceptionHandler() { }
    virtual ~GExceptionHandler() throw() { }
    virtual const char* what() const throw();
protected:
    std::string m_origin;
    std::string m_message;
};


/***********************************************************************//**
 * @class GException
 *
 * @brief Interface for exceptions
 *
 * The exception class is the master class that is thrown in case of
 * exceptions.
 ***************************************************************************/
class GException : public GExceptionHandler {
public:

    // --- LOGIC EXCEPTIONS (what client could have tested) ---

    // Invalid value
    class invalid_value : public GExceptionHandler {
    public:
        invalid_value(const std::string& origin,
                      const std::string& message);
    };

    // Invalid argument
    class invalid_argument : public GExceptionHandler {
    public:
        invalid_argument(const std::string& origin,
                         const std::string& message);
        invalid_argument(const std::string& origin,
                         const std::string& argument,
                         const std::string& message);
    };

    // Invalid return value
    class invalid_return_value : public GExceptionHandler {
    public:
        invalid_return_value(const std::string& origin,
                             const std::string& message);
    };

    // Out of range
    class out_of_range : public GExceptionHandler {
    public:
        out_of_range(const std::string& origin,
                     const std::string& what,
                     const int&         index,
                     const int&         elements,
                     const std::string& message = "");
    };

    // FITS error
    class fits_error : public GExceptionHandler {
    public:
        fits_error(const std::string& origin,
		           const int&         status,
				   const std::string& message = "");
    };


    // --- RUNTIME EXCEPTIONS (not testable by client) ---

    // Runtime error
    class runtime_error : public GExceptionHandler {
    public:
        runtime_error(const std::string& origin,
                      const std::string& message = "");
    };

    //underflow_error

    //overflow_error

    // File error
    class file_error : public GExceptionHandler {
    public:
        file_error(const std::string& origin,
                   const std::string& message = "");
    };

    // Feature not implemented
    class feature_not_implemented : public GExceptionHandler {
    public:
        feature_not_implemented(const std::string& origin,
                                const std::string& message = "");
    };


    // --- TEST EXCEPTIONS (only used by test suite) ---

    // Signal nested try
    class test_nested_try_error : public GExceptionHandler {
        public:
            test_nested_try_error(const std::string& origin,
                                  const std::string& message = "");
    };

    // Signal test failure
    class test_failure : public GExceptionHandler {
        public:
            test_failure(const std::string& origin,
                         const std::string& message = "");
    };

    // Signal test error
    class test_error : public GExceptionHandler {
        public:
            test_error(const std::string& origin,
                       const std::string& message = "");
    };

};

#endif /* GEXCEPTION_HPP */
