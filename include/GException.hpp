/***************************************************************************
 *                   GException.hpp  -  exception handler                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2006-2011 by Jurgen Knodlseder                           *
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
 * @brief Exception handler interface definiton.
 * @author J. Knodlseder
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
 * @brief Interface for exceptions.
 *
 * This is the object that it thrown in case of an exception.
 ***************************************************************************/
class GException : public GExceptionHandler {
public:

    // Feature not implemented
    class feature_not_implemented : public GExceptionHandler {
    public:
        feature_not_implemented(std::string origin,
                                std::string message = "");
    };

    // Invalid argument
    class invalid_argument : public GExceptionHandler {
    public:
        invalid_argument(std::string origin,
                         std::string message = "");
    };

    // Bad type
    class bad_type : public GExceptionHandler {
    public:
        bad_type(std::string origin, std::string message = "");
    };

    // Memory allocation exception class
    class mem_alloc : public GExceptionHandler {
    public:
        mem_alloc(std::string origin, unsigned num);
    };

    // Not enough nodes
    class not_enough_nodes : public GExceptionHandler {
    public:
        not_enough_nodes(std::string origin, int num);
    };

    // File not found
    class file_not_found : public GExceptionHandler {
    public:
        file_not_found(std::string origin, std::string filename);
    };

    // File open error
    class file_open_error : public GExceptionHandler {
    public:
        file_open_error(std::string origin, std::string filename);
    };

    // Empty object exception class
    class empty : public GExceptionHandler {
    public:
        empty(std::string origin);
    };

    // Out of range
    class out_of_range : public GExceptionHandler {
    public:
        out_of_range(std::string origin, int inx, int min, int max);
        out_of_range(std::string origin, double value, double min, double max);
        out_of_range(std::string origin, int inx, int elements);
        out_of_range(std::string origin, int row, int col, int rows, int cols);
    };

    // Vector - Vector mismatch
    class vector_mismatch : public GExceptionHandler {
    public:
        vector_mismatch(std::string origin, int size1, int size2);
    };

    // Cross product only defined for 3-element vectors
    class vector_bad_cross_dim : public GExceptionHandler {
    public:
        vector_bad_cross_dim(std::string origin, int elements);
    };

    // Vector - Matrix mismatch
    class matrix_vector_mismatch : public GExceptionHandler {
    public:
        matrix_vector_mismatch(std::string origin, int num, int rows, int cols);
    };

    // Matrix dimensions mismatch
    class matrix_mismatch : public GExceptionHandler {
    public:
        matrix_mismatch(std::string origin, int rows1, int cols1, int rows2, 
                        int cols2);
    };

    // Matrix not square
    class matrix_not_square : public GExceptionHandler {
    public:
        matrix_not_square(std::string origin, int rows, int cols);
    };

    // Matrix not positive definite
    class matrix_not_pos_definite : public GExceptionHandler {
    public:
        matrix_not_pos_definite(std::string origin, int row, double sum);
    };

    // Matrix not symmetric
    class matrix_not_symmetric : public GExceptionHandler {
    public:
        matrix_not_symmetric(std::string origin, int cols, int rows);
    };

    // Matrix not factorised
    class matrix_not_factorised : public GExceptionHandler {
    public:
        matrix_not_factorised(std::string origin, std::string type);
    };

    // All matrix elements are zero
    class matrix_zero : public GExceptionHandler  {
    public:
        matrix_zero(std::string origin);
    };

    // Invalid ordering scheme
    class invalid_order : public GExceptionHandler {
    public:
        invalid_order(std::string origin, int order, int min_order, 
                      int max_order);
    };


    // FITS exceptions
    class fits_error : public GExceptionHandler {
    public:
        fits_error(std::string origin,
		           int         status,
				   std::string message = "");
    };
    class fits_open_error : public GExceptionHandler {
    public:
        fits_open_error(std::string origin,
		                std::string filename,
		                int         status,
						std::string message = "");
    };
    class fits_file_exist : public GExceptionHandler {
    public:
        fits_file_exist(std::string origin,
		                std::string filename,
						int         status = 0);
    };
    class fits_file_not_open : public GExceptionHandler {
    public:
        fits_file_not_open(std::string origin,
		                   std::string filename);
    };
    class fits_already_opened : public GExceptionHandler {
    public:
        fits_already_opened(std::string origin,
		                    std::string filename);
    };
    class fits_key_not_found : public GExceptionHandler {
    public:
        fits_key_not_found(std::string origin,
		                   std::string keyname, 
                           int         status = 0);
    };
    class fits_column_not_found : public GExceptionHandler {
    public:
        fits_column_not_found(std::string origin,
		                      std::string colname, 
                              int         status = 0);
    };
    class fits_no_header : public GExceptionHandler {
    public:
        fits_no_header(std::string origin,
		               std::string message,
					   int         status = 0);
    };
    class fits_no_data : public GExceptionHandler {
    public:
        fits_no_data(std::string origin,
		             std::string message,
					 int         status = 0);
    };
    class fits_hdu_not_found : public GExceptionHandler {
    public:
        fits_hdu_not_found(std::string origin,
		                   std::string extname,
                           int         status = 0);
        fits_hdu_not_found(std::string origin,
		                   int         extno,
						   int         status = 0);
    };
    class fits_hdu_not_image : public GExceptionHandler {
    public:
        fits_hdu_not_image(std::string origin,
		                   std::string extname,
						   int         type);
    };
    class fits_hdu_not_table : public GExceptionHandler {
    public:
        fits_hdu_not_table(std::string origin,
		                   std::string extname,
						   int         type);
    };
    class fits_unknown_HDU_type : public GExceptionHandler {
    public:
        fits_unknown_HDU_type(std::string origin,
		                      int         type);
    };
    class fits_invalid_type : public GExceptionHandler {
    public:
        fits_invalid_type(std::string origin,
		                  std::string message);
    };
    class fits_unknown_tabtype : public GExceptionHandler {
    public:
        fits_unknown_tabtype(std::string origin,
		                     int         type);
    };
    class fits_unknown_coltype : public GExceptionHandler {
    public:
        fits_unknown_coltype(std::string origin,
		                     std::string colname,
							 int         type);
    };
    class fits_bad_col_length : public GExceptionHandler {
    public:
        fits_bad_col_length(std::string origin,
		                    int         length,
							int         rows);
    };
    class fits_bad_bitpix : public GExceptionHandler {
    public:
        fits_bad_bitpix(std::string origin,
		                int         bitpix);
    };
    class fits_wrong_image_operator : public GExceptionHandler {
    public:
        fits_wrong_image_operator(std::string origin,
		                          int         naxis,
								  int         nargs);
    };
    class fits_invalid_row : public GExceptionHandler {
    public:
        fits_invalid_row(std::string origin,
                         int         row,
                         int         nrows,
                         std::string message = "");
    };
    class fits_invalid_nrows : public GExceptionHandler {
    public:
        fits_invalid_nrows(std::string origin,
                           int         nrows,
                           int         max_rows,
                           std::string message = "");
    };
    class fits_inconsistent_tdim : public GExceptionHandler {
    public:
        fits_inconsistent_tdim(std::string      origin,
                               std::vector<int> tdim,
                               int              number,
                               std::string      message = "");
    };


    // GSkymap exceptions
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


    // Application exceptions
    class app_error : public GExceptionHandler {
    public:
        app_error(std::string origin,
                  std::string message = "");
    };
    class par_file_not_found : public GExceptionHandler {
    public:
        par_file_not_found(std::string origin,
                           std::string filename,
                           std::string message = "");
    };
    class par_file_open_error : public GExceptionHandler {
    public:
        par_file_open_error(std::string origin,
                            std::string filename,
                            std::string message = "");
    };
    class home_not_found : public GExceptionHandler {
    public:
        home_not_found(std::string origin,
                       std::string message = "");
    };
    class could_not_create_pfiles : public GExceptionHandler {
    public:
        could_not_create_pfiles(std::string origin,
                                std::string home,
                                std::string message = "");
    };
    class pfiles_not_accessible : public GExceptionHandler {
    public:
        pfiles_not_accessible(std::string origin,
                              std::string home,
                              std::string message = "");
    };
    class par_file_syntax_error : public GExceptionHandler {
    public:
        par_file_syntax_error(std::string origin,
                              std::string home,
                              std::string message = "");
    };
    class par_error : public GExceptionHandler {
    public:
        par_error(std::string origin,
                  std::string name,
                  std::string message = "");
    };
    class bad_cmdline_argument : public GExceptionHandler {
    public:
        bad_cmdline_argument(std::string origin,
                             std::string arg,
                             std::string message = "");
    };


    // Observation exceptions
    class no_response : public GExceptionHandler {
    public:
        no_response(std::string origin, std::string message = "");
    };
    class no_roi : public GExceptionHandler {
    public:
        no_roi(std::string origin, std::string message = "");
    };
    class no_events : public GExceptionHandler {
    public:
        no_events(std::string origin, std::string message = "");
    };
    class no_list : public GExceptionHandler {
    public:
        no_list(std::string origin, std::string message = "");
    };
    class no_cube : public GExceptionHandler {
    public:
        no_cube(std::string origin, std::string message = "");
    };
    class gradient_par_mismatch : public GExceptionHandler {
    public:
        gradient_par_mismatch(std::string origin, int nsize, int npars);
    };
    class caldb_not_found : public GExceptionHandler {
    public:
        caldb_not_found(std::string origin, std::string home, 
                        std::string message = "");
    };
    class rsp_invalid_type : public GExceptionHandler {
    public:
        rsp_invalid_type(std::string origin, std::string type);
    };
    class gti_invalid : public GExceptionHandler {
    public:
        gti_invalid(std::string origin, GTime tstart, GTime tstop,
                    std::string message = "");
    };
    class erange_invalid : public GExceptionHandler {
    public:
        erange_invalid(std::string origin, double emin, 
                       double emax, std::string message = "");
    };
    class invalid_statistics : public GExceptionHandler {
    public:
        invalid_statistics(std::string origin, std::string statistics,
                           std::string message = "");
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


    // Model exceptions
    class model_invalid : public GExceptionHandler {
    public:
        model_invalid(std::string origin, std::string type,
                      std::string message = "");
    };
    class model_invalid_spatial : public GExceptionHandler {
    public:
        model_invalid_spatial(std::string origin, std::string type,
                              std::string message = "");
    };
    class model_invalid_spectral : public GExceptionHandler {
    public:
        model_invalid_spectral(std::string origin, std::string type,
                               std::string message = "");
    };
    class model_invalid_temporal : public GExceptionHandler {
    public:
        model_invalid_temporal(std::string origin, std::string type,
                               std::string message = "");
    };
    class model_invalid_parnum : public GExceptionHandler {
    public:
        model_invalid_parnum(std::string origin, GXmlElement xml,
                             std::string message = "");
    };
    class model_invalid_parnames : public GExceptionHandler {
    public:
        model_invalid_parnames(std::string origin, GXmlElement xml,
                               std::string message = "");
    };
    class model_invalid_parscale : public GExceptionHandler {
    public:
        model_invalid_parscale(std::string origin, double scale,
                               std::string message = "");
        model_invalid_parscale(std::string origin, GXmlElement xml,
                               std::string message = "");
    };
    class file_function_data : public GExceptionHandler {
    public:
        file_function_data(std::string origin, std::string filename,
                           int num, std::string message = "");
    };
    class file_function_columns : public GExceptionHandler {
    public:
        file_function_columns(std::string origin, std::string filename,
                              int num, std::string message = "");
    };
    class file_function_value : public GExceptionHandler {
    public:
        file_function_value(std::string origin, std::string filename,
                            double value, std::string message = "");
    };
    class par_not_found : public GExceptionHandler {
    public:
        par_not_found(std::string origin, std::string name,
                      std::string message = "");
    };
    class model_not_found : public GExceptionHandler {
    public:
        model_not_found(std::string origin, std::string name,
                        std::string message = "");
    };
    class no_point_source : public GExceptionHandler {
    public:
        no_point_source(std::string origin, std::string name,
                        std::string message = "");
    };
    class no_extended_source : public GExceptionHandler {
    public:
        no_extended_source(std::string origin, std::string name,
                           std::string message = "");
    };

    // CSV exceptions
    class csv_bad_columns : public GExceptionHandler {
    public:
        csv_bad_columns(std::string origin, std::string filename,
                        int rows, int cols, int elements,
                        std::string message = "");
    };

};

#endif /* GEXCEPTION_HPP */
