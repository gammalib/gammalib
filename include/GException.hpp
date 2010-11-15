/***************************************************************************
 *                   GException.hpp  -  exception handler                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2006-2010 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
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

    // Memory allocation exception class
    class mem_alloc : public GExceptionHandler {
    public:
        mem_alloc(std::string origin, unsigned num);
    };
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

    // General FITS error
    class fits_error : public GExceptionHandler {
    public:
        fits_error(std::string origin, int status, std::string message = "");
    };

    // FITS file open error
    class fits_open_error : public GExceptionHandler {
    public:
        fits_open_error(std::string origin, std::string filename, int status);
    };

    // FITS file exists already
    class fits_file_exist : public GExceptionHandler {
    public:
        fits_file_exist(std::string origin, std::string filename, int status);
    };

    // FITS file has already been opened
    class fits_already_opened : public GExceptionHandler {
    public:
        fits_already_opened(std::string origin, std::string filename);
    };

    // FITS keyword not found error
    class fits_key_not_found : public GExceptionHandler {
    public:
        fits_key_not_found(std::string origin, std::string keyname, 
                           int status = 0);
    };

    // FITS column not found error
    class fits_column_not_found : public GExceptionHandler {
    public:
        fits_column_not_found(std::string origin, std::string colname, 
                              int status = 0);
    };

    // FITS has no header
    class fits_no_header : public GExceptionHandler {
    public:
        fits_no_header(std::string origin, std::string message, int status = 0);
    };

    // FITS has no data
    class fits_no_data : public GExceptionHandler {
    public:
        fits_no_data(std::string origin, std::string message, int status = 0);
    };

    // FITS HDU not found error
    class fits_hdu_not_found : public GExceptionHandler {
    public:
        fits_hdu_not_found(std::string origin, std::string extname, 
                           int status = 0);
        fits_hdu_not_found(std::string origin, int extno, int status = 0);
    };

    // FITS unknown HDU type
    class fits_unknown_HDU_type : public GExceptionHandler {
    public:
        fits_unknown_HDU_type(std::string origin, int type);
    };

    // FITS HDU is not a table
    class fits_HDU_not_a_table : public GExceptionHandler {
    public:
        fits_HDU_not_a_table(std::string origin, int type);
    };

    // FITS invalid type
    class fits_invalid_type : public GExceptionHandler {
    public:
        fits_invalid_type(std::string origin, std::string message);
    };

    // FITS unknown table type
    class fits_unknown_tabtype : public GExceptionHandler {
    public:
        fits_unknown_tabtype(std::string origin, int type);
    };

    // FITS unknown column type
    class fits_unknown_coltype : public GExceptionHandler {
    public:
        fits_unknown_coltype(std::string origin, std::string colname, int type);
    };

    // FITS bad column length
    class fits_bad_col_length : public GExceptionHandler {
    public:
        fits_bad_col_length(std::string origin, int length, int rows);
    };

    // FITS bad bitpix value
    class fits_bad_bitpix : public GExceptionHandler {
    public:
        fits_bad_bitpix(std::string origin, int bitpix);
    };

    // FITS bad image operator
    class fits_wrong_image_operator : public GExceptionHandler {
    public:
        fits_wrong_image_operator(std::string origin, int naxis, int nargs);
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


    // Application exceptions
    class par_file_not_found : public GExceptionHandler {
    public:
        par_file_not_found(std::string origin, std::string filename,
                           std::string message = "");
    };
    class par_file_open_error : public GExceptionHandler {
    public:
        par_file_open_error(std::string origin, std::string filename,
                            std::string message = "");
    };
    class home_not_found : public GExceptionHandler {
    public:
        home_not_found(std::string origin, std::string message = "");
    };
    class could_not_create_pfiles : public GExceptionHandler {
    public:
        could_not_create_pfiles(std::string origin, std::string home,
                                std::string message = "");
    };
    class pfiles_not_accessible : public GExceptionHandler {
    public:
        pfiles_not_accessible(std::string origin, std::string home,
                              std::string message = "");
    };
    class par_file_syntax_error : public GExceptionHandler {
    public:
        par_file_syntax_error(std::string origin, std::string home,
                              std::string message = "");
    };
    class par_error : public GExceptionHandler {
    public:
        par_error(std::string origin, std::string message = "");
    };
    class bad_cmdline_argument : public GExceptionHandler {
    public:
        bad_cmdline_argument(std::string origin, std::string arg,
                             std::string message = "");
    };


    // Observation exceptions
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
    class roi_invalid : public GExceptionHandler {
    public:
        roi_invalid(std::string origin, std::string message = "");
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


    // GModel exceptions
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

};

#endif /* GEXCEPTION_HPP */
