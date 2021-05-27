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
#include <vector>                             // string
#include <sstream>                            // ostringstream
#include <stdexcept>                          // exception
#include "GTime.hpp"
#include "GXmlElement.hpp"

/* __ Prototypes ________________________________________________________ */
namespace gammalib {
    void check_energy_interval(const std::string& origin,
                               const GEnergy&     emin,
                               const GEnergy&     emax);
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




    // --- OLD EXCEPTIONS ---

    // GSkyMap exceptions
    class skymap : public GExceptionHandler {
    public:
        skymap(std::string origin, std::string message = "");
    };
    class skymap_bad_par : public GExceptionHandler {
    public:
        skymap_bad_par(std::string origin, int par, std::string message = "");
    };
    class skymap_bad_size : public GExceptionHandler {
    public:
        skymap_bad_size(std::string origin, int size, int expected,
                        std::string message = "");
    };
    class skymap_bad_ctype : public GExceptionHandler {
    public:
        skymap_bad_ctype(std::string origin, std::string ctype1,
                         std::string ctype2, std::string message = "");
    };
    class skymap_bad_image_dim : public GExceptionHandler {
    public:
        skymap_bad_image_dim(std::string origin, int naxis,
                             std::string message = "");
    };


    // GWcs exceptions
    class wcs : public GExceptionHandler {
    public:
        wcs(std::string origin, std::string message = "");
    };
    class wcs_invalid : public GExceptionHandler {
    public:
        wcs_invalid(std::string origin, std::string wcs, 
                    std::string message = "");
    };
    class wcs_bad_coords : public GExceptionHandler {
    public:
        wcs_bad_coords(std::string origin, std::string coordsys);
    };
    class wcs_no_proj_fct : public GExceptionHandler {
    public:
        wcs_no_proj_fct(std::string origin, std::string message = "");
    };
    class wcs_hpx_bad_nside : public GExceptionHandler {
    public:
        wcs_hpx_bad_nside(std::string origin, int nside);
    };
    class wcs_hpx_bad_ordering : public GExceptionHandler {
    public:
        wcs_hpx_bad_ordering(std::string origin, std::string ordering);
    };
    class wcs_singular_matrix : public GExceptionHandler {
    public:
        wcs_singular_matrix(std::string origin, int naxis,
                            const std::vector<double>& mat);
    };
    class wcs_invalid_parameter : public GExceptionHandler {
    public:
        wcs_invalid_parameter(std::string origin, std::string message = "");
    };
    class wcs_invalid_x_y : public GExceptionHandler {
    public:
        wcs_invalid_x_y(std::string origin, int num, std::string message = "");
    };
    class wcs_invalid_phi_theta : public GExceptionHandler {
    public:
        wcs_invalid_phi_theta(std::string origin, int num, std::string message = "");
    };


    // XML exceptions
    class xml_syntax_error : public GExceptionHandler {
    public:
        xml_syntax_error(std::string origin, std::string segment,
                         std::string message = "");
    };
    class xml_attribute_value : public GExceptionHandler {
    public:
        xml_attribute_value(std::string origin, std::string value);
    };
    class xml_bad_node_type : public GExceptionHandler {
    public:
        xml_bad_node_type(std::string origin, std::string type,
                          std::string message = "");
    };
    class xml_name_not_found : public GExceptionHandler {
    public:
        xml_name_not_found(std::string origin, std::string type,
                           std::string message = "");
    };
    class xml_invalid_parnum : public GExceptionHandler {
    public:
        xml_invalid_parnum(std::string origin, GXmlElement xml,
                           std::string message = "");
    };
    class xml_invalid_parnames : public GExceptionHandler {
    public:
        xml_invalid_parnames(std::string origin, GXmlElement xml,
                             std::string message = "");
    };

    // CSV exceptions
    class csv_bad_columns : public GExceptionHandler {
    public:
        csv_bad_columns(std::string origin, std::string filename,
                        int rows, int cols, int elements,
                        std::string message = "");
    };

    // Test exceptions
    class test_nested_try_error : public GExceptionHandler {
        public:
            test_nested_try_error(std::string origin, std::string message = "");
    };
    
    class test_failure : public GExceptionHandler {
        public:
            test_failure(std::string origin, std::string message = "");
    };
    
    class test_error : public GExceptionHandler {
        public:
            test_error(std::string origin, std::string message = "");
    };

};

#endif /* GEXCEPTION_HPP */
