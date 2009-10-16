/***************************************************************************
 *                   GException.hpp  -  exception handler                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2006-2009 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef GEXCEPTION_HPP
#define GEXCEPTION_HPP

/* __ Includes ___________________________________________________________ */
#include <string>                             // string
#include <sstream>                            // ostringstream
#include <stdexcept>                          // exception


/***************************************************************************
 *                  Exception handler base class definition                *
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


/***************************************************************************
 *                         Exception class definition                      *
 ***************************************************************************/
class GException : public GExceptionHandler {
public:

    // Memory allocation exception class
    class mem_alloc : public GExceptionHandler {
    public:
        mem_alloc(std::string origin, unsigned num);
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

    // Matrix not rectangular
    class matrix_not_rectangular : public GExceptionHandler {
    public:
        matrix_not_rectangular(std::string origin, int rows, int cols);
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
        invalid_order(std::string origin, int order, int min_order, int max_order);
    };

    // General FITS error
    class fits_error : public GExceptionHandler {
    public:
        fits_error(std::string origin, int status);
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
        fits_key_not_found(std::string origin, std::string keyname, int status = 0);
    };

    // FITS column not found error
    class fits_column_not_found : public GExceptionHandler {
    public:
        fits_column_not_found(std::string origin, std::string colname, int status = 0);
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
        fits_hdu_not_found(std::string origin, std::string extname, int status = 0);
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
        fits_unknown_coltype(std::string origin, int type);
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

    // Response invalid response type
    class rsp_invalid_type : public GExceptionHandler {
    public:
        rsp_invalid_type(std::string origin, std::string type);
    };

    // General Healpix error
    class healpix : public GExceptionHandler {
    public:
        healpix(std::string origin, std::string message);
    };

    // Healpix resolution invalid
    class healpix_bad_nside : public GExceptionHandler {
    public:
        healpix_bad_nside(std::string origin, int nside);
    };

    // Healpix scheme invalid
    class healpix_bad_scheme : public GExceptionHandler {
    public:
        healpix_bad_scheme(std::string origin, std::string scheme);
    };

    // Healpix coordinate system invalid
    class healpix_bad_coords : public GExceptionHandler {
    public:
        healpix_bad_coords(std::string origin, std::string coordsys);
    };

    // Invalid object release
    class gradient_par_mismatch : public GExceptionHandler {
    public:
        gradient_par_mismatch(std::string origin, int nsize, int npars);
    };

};

#endif /* GEXCEPTION_HPP */
