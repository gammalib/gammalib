/***************************************************************************
 *                   GException.cpp  -  exception handler                  *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2006 by Jurgen Knodlseder                   *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

/* __ Includes ___________________________________________________________ */
#include "GException.hpp"
#include "GFitsCfitsio.hpp"

/* __ Namespaces _________________________________________________________ */
using namespace std;


/***************************************************************************
 *                       Private conversion functions                      *
 ***************************************************************************/
std::string str(int value)
{
    ostringstream s_value;
    s_value << value;
    return  s_value.str();
}
std::string str(double value)
{
    ostringstream s_value;
    s_value << scientific << value;
    return  s_value.str();
}


/***************************************************************************
 *                  Exception handler base class definition                *
 ***************************************************************************/
const char* GExceptionHandler::what() const throw()
{
    string message = "*** ERROR in " + m_origin + ": " + m_message;
    return message.c_str();
}


/***************************************************************************
 *                        Memory allocation exception                      *
 ***************************************************************************/
GException::mem_alloc::mem_alloc(string origin, unsigned num)
{
    m_origin  = origin;
    m_message = "Memory allocation error (" + str((int)num) + " elements)";
}


/***************************************************************************
 *                           Empty object exception                        *
 ***************************************************************************/
GException::empty::empty(string origin)
{
    m_origin  = origin;
    m_message = "Zero-size allocation";
}


/***************************************************************************
 *                            Index out of range                           *
 ***************************************************************************/
GException::out_of_range::out_of_range(string origin, int inx, int min, int max)
{
    m_origin  = origin;
    m_message = "Index (" + str(inx) + ") out of range [" + str(min) +
                "," + str(max) + "]";
}


/***************************************************************************
 *                          Vector index out of range                      *
 ***************************************************************************/
GException::out_of_range::out_of_range(string origin, int inx, int elements)
{
    m_origin = origin;
    if (elements > 0) {
        m_message = "Vector index (" + str(inx) + ") out of range [0," +
                    str(elements-1) + "]";
    }
    else {
        m_message = "Empty vector";
    }
}


/***************************************************************************
 *                      Matrix row or column out of range                  *
 ***************************************************************************/
GException::out_of_range::out_of_range(string origin, int row, int col, int rows, int cols)
{
    m_origin  = origin;
    m_message = "Matrix element (" + str(row) + "," + str(col) +
                ") out of range ([0," + str(rows) + "], [0," +
                str(cols) + "])";
}


/***************************************************************************
 *                          Vector dimensions differ                       *
 ***************************************************************************/
GException::vector_mismatch::vector_mismatch(string origin, int size1, int size2)
{
    m_origin  = origin;
    m_message = "Vector dimensions differ (" + str(size1) + " <-> " +
                str(size2) + ")";
}


/***************************************************************************
 *                   Invalid vector dimension for cross product            *
 ***************************************************************************/
GException::vector_bad_cross_dim::vector_bad_cross_dim(string origin, int elements)
{
    m_origin  = origin;
    m_message = "Vector cross product only defined for 3 dimensions but vector size is " + 
                str(elements); 
}


/***************************************************************************
 *                     Mismatch between matrix and vector                  *
 ***************************************************************************/
GException::matrix_vector_mismatch::matrix_vector_mismatch(string origin, int num, int rows, int cols)
{
    m_origin  = origin;
    m_message = "Vector dimension [" + str(num) + 
                "] is incompatible with matrix size [" + 
                str(rows) + "," + str(cols) + "]";
}


/***************************************************************************
 *                       Matrix mismatch in operation                      *
 ***************************************************************************/
GException::matrix_mismatch::matrix_mismatch(string origin, int rows1, int cols1, int rows2, int cols2)
{
    m_origin  = origin;
    m_message = "Matrix mismatch: M1(" + str(rows1) + "," + str(cols1) +
                ") incompatible with M2(" + str(rows2) + "," + str(cols2) + ")";
}


/***************************************************************************
 *                           Matrix not rectangular                        *
 ***************************************************************************/
GException::matrix_not_rectangular::matrix_not_rectangular(string origin, int rows, int cols)
{
    m_origin  = origin;
    m_message = "Matrix is not rectangular [" + str(rows) + "," + str(cols) + "]";
}


/***************************************************************************
 *                     Matrix is not positive definite                     *
 ***************************************************************************/
GException::matrix_not_pos_definite::matrix_not_pos_definite(string origin, int row, double sum)
{
    m_origin  = origin;
    m_message = "Matrix is not positive definite (sum " + str(sum) + 
                " occured in row/column " + str(row) + ")";
}


/***************************************************************************
 *                         Matrix is not symmetric                         *
 ***************************************************************************/
GException::matrix_not_symmetric::matrix_not_symmetric(string origin, int cols, int rows)
{
    m_origin  = origin;
    m_message = "Matrix is not symmetric [" + str(rows) + "," + str(cols) + "]";
}


/***************************************************************************
 *                       Matrix has not been factorised                    *
 ***************************************************************************/
GException::matrix_not_factorised::matrix_not_factorised(string origin, string type)
{
    m_origin  = origin;
    m_message = "Matrix has not been factorised using " + type;
}


/***************************************************************************
 *                       All matrix elements are zero                      *
 ***************************************************************************/
GException::matrix_zero::matrix_zero(string origin)
{
    m_origin = origin;
    m_message = "All matrix elements are zero";
}


/***************************************************************************
 *                     Invalid ordering scheme requested                   *
 ***************************************************************************/
GException::invalid_order::invalid_order(string origin, int order, int min_order, int max_order)
{
    m_origin  = origin;
    m_message = "Invalid ordering type " + str(order) + 
                "requested; must be comprised in [" + str(min_order) +
                "," + str(max_order) + "]";
}


/***********************************************************************//**
 * @brief General FITS error
 *
 * @param[in] origin Method that throws the error.
 * @param[in] status cfitsio status.
 ***************************************************************************/
GException::fits_error::fits_error(string origin, int status)
{
    m_origin  = origin;
    ostringstream s_error;
    char err_text[31];
    __ffgerr(status, err_text);
    s_error << err_text << " (status=" << status << ")";
    m_message = s_error.str();
}


/***********************************************************************//**
 * @brief FITS error: unable to open FITS file
 *
 * @param[in] origin Method that throws the error.
 * @param[in] filename Name of the file for which opening was attempted.
 * @param[in] status cfitsio status.
 ***************************************************************************/
GException::fits_open_error::fits_open_error(string origin,
                                             string filename,
                                             int    status)
{
    m_origin  = origin;
    m_message = "Unable to open FITS file '" + filename + "'";
    if (status != 0)
        m_message += " (status=" + str(status) + ")";
}


/***********************************************************************//**
 * @brief FITS error: attempted to overwrite FITS file
 *
 * @param[in] origin Method that throws the error.
 * @param[in] filename Name of the file for which overwrite attempt was made.
 * @param[in] status cfitsio status.
 ***************************************************************************/
GException::fits_file_exist::fits_file_exist(string origin,
                                             string filename,
                                             int    status)
{
    m_origin  = origin;
    m_message = "Attempted to overwrite FITS file '" + filename + "'";
    if (status != 0)
        m_message += " (status=" + str(status) + ")";
}


/***********************************************************************//**
 * @brief FITS error: file already open
 *
 * @param[in] origin Method that throws the error.
 * @param[in] filename Name of file that is already open.
 ***************************************************************************/
GException::fits_already_opened::fits_already_opened(string origin,
                                                     string filename)
{
    m_origin  = origin;
    m_message = "FITS file '" + filename + "' is already open";
}


/***********************************************************************//**
 * @brief FITS error: Keyword not in header
 *
 * @param[in] origin Method that throws the error.
 * @param[in] keyname Name of keyword that was not found.
 * @param[in] status cfitsio status.
 ***************************************************************************/
GException::fits_key_not_found::fits_key_not_found(string origin,
                                                   string keyname,
                                                   int    status)
{
    m_origin  = origin;
    m_message = "Keyword '" + keyname + "' not found in header";
    if (status != 0)
        m_message += " (status=" + str(status) + ")";
}


/***********************************************************************//**
 * @brief FITS error: Table column not found
 *
 * @param[in] origin Method that throws the error.
 * @param[in] colname Name of the table column that was not found.
 * @param[in] status cfitsio status.
 ***************************************************************************/
GException::fits_column_not_found::fits_column_not_found(string origin,
                                                         string colname,
                                                         int    status)
{
    m_origin  = origin;
    m_message = "Column '" + colname + "' not found in table";
    if (status != 0)
        m_message += " (status=" + str(status) + ")";
}


/***********************************************************************//**
 * @brief FITS error: No header
 *
 * @param[in] origin Method that throws the error.
 * @param[in] message Error message.
 * @param[in] status cfitsio status.
 ***************************************************************************/
GException::fits_no_header::fits_no_header(string origin,
                                           string message,
                                           int    status)
{
    m_origin  = origin;
    m_message = message;
    if (status != 0)
        m_message += " (status=" + str(status) + ")";
}


/***********************************************************************//**
 * @brief FITS error: No data
 *
 * @param[in] origin Method that throws the error.
 * @param[in] message Error message.
 * @param[in] status cfitsio status.
 ***************************************************************************/
GException::fits_no_data::fits_no_data(string origin,
                                       string message,
                                       int    status)
{
    m_origin  = origin;
    m_message = message;
    if (status != 0)
        m_message += " (status=" + str(status) + ")";
}


/***********************************************************************//**
 * @brief FITS error: HDU not found in FITS file
 *
 * @param[in] origin Method that throws the error.
 * @param[in] extname Name of the extension that was not found.
 * @param[in] status cfitsio status.
 ***************************************************************************/
GException::fits_hdu_not_found::fits_hdu_not_found(string origin,
                                                   string extname,
                                                   int    status)
{
    m_origin  = origin;
    m_message = "HDU '" + extname + "' not found in FITS file";
    if (status != 0)
        m_message += " (status=" + str(status) + ")";
}


/***********************************************************************//**
 * @brief FITS error: HDU not found in FITS file
 *
 * @param[in] origin Method that throws the error.
 * @param[in] extno Number of the extension that was not found.
 * @param[in] status cfitsio status.
 ***************************************************************************/
GException::fits_hdu_not_found::fits_hdu_not_found(string origin,
                                                   int    extno,
                                                   int    status)
{
    m_origin  = origin;
    m_message = "HDU number '" + str(extno) + "' not found in FITS file";
    if (status != 0)
        m_message += " (status=" + str(status) + ")";
}


/***********************************************************************//**
 * @brief FITS error: HDU of unknown type found
 *
 * @param[in] origin Method that throws the error.
 * @param[in] type Specified HDU type.
 ***************************************************************************/
GException::fits_unknown_HDU_type::fits_unknown_HDU_type(string origin,
                                                         int    type)
{
    m_origin  = origin;
    m_message = "HDU type '" + str(type) + "' is not known";
}


/***********************************************************************//**
 * @brief FITS error: HDU is not a table
 *
 * @param[in] origin Method that throws the error.
 * @param[in] type Specified HDU type.
 ***************************************************************************/
GException::fits_HDU_not_a_table::fits_HDU_not_a_table(string origin,
                                                       int    type)
{
    m_origin  = origin;
    m_message = "HDU is not of type 'table' (type=" + str(type) + ")";
}


/***********************************************************************//**
 * @brief FITS error: Invalid type
 *
 * @param[in] origin Method that throws the error.
 * @param[in] message Error message.
 ***************************************************************************/
GException::fits_invalid_type::fits_invalid_type(string origin,
                                                 string message)
{
    m_origin  = origin;
    m_message = message;
}


/***********************************************************************//**
 * @brief FITS error: Table type is unknow
 *
 * @param[in] origin Method that throws the error.
 * @param[in] type Specified table type.
 ***************************************************************************/
GException::fits_unknown_tabtype::fits_unknown_tabtype(string origin,
                                                       int    type)
{
    m_origin  = origin;
    m_message = "Table type '" + str(type) + "' is unknown";
}


/***********************************************************************//**
 * @brief FITS error: Column type is unknown
 *
 * @param[in] origin Method that throws the error.
 * @param[in] type Specified column type.
 ***************************************************************************/
GException::fits_unknown_coltype::fits_unknown_coltype(string origin,
                                                       int    type)
{
    m_origin  = origin;
    m_message = "Column type '" + str(type) + "' is unknown";
}


/***********************************************************************//**
 * @brief FITS error: Bad column length
 *
 * @param[in] origin Method that throws the error.
 * @param[in] length Length of the column.
 * @param[in] rows Number of rows in table.
 ***************************************************************************/
GException::fits_bad_col_length::fits_bad_col_length(string origin,
                                                     int    length,
                                                     int    rows)
{
    m_origin  = origin;
    m_message = "Column length '" + str(length) + "' is not compatible"
                " with the number of rows '" + str(rows) + "' in the"
                " table";
}


/***********************************************************************//**
 * @brief FITS error: invalid number of bits per pixel
 *
 * @param[in] origin Method that throws the error.
 * @param[in] bitpix Bitpix value that was not 8,16,32,64,-32, or -64.
 ***************************************************************************/
GException::fits_bad_bitpix::fits_bad_bitpix(string origin,
                                             int    bitpix)
{
    m_origin  = origin;
    m_message = "Invalid number of bits per pixel (bitpix="+str(bitpix)+")";
}


/***********************************************************************//**
 * @brief FITS error: wrong image operator has been used
 *
 * @param[in] origin Method that throws the error.
 * @param[in] naxis Dimension of image.
 * @param[in] nargs Number of arguments of the image operator.
 ***************************************************************************/
GException::fits_wrong_image_operator::fits_wrong_image_operator(string origin,
                                                                 int    naxis,
                                                                 int    nargs)
{
    m_origin  = origin;
    m_message = "Wrong image pixel access operator has been used (dimension=" +
                str(naxis) + " <=> arguments=" + str(nargs) + ")";
}


/***********************************************************************//**
 * @brief Response error: invalid response type specified
 *
 * @param[in] origin Method that throws the error.
 * @param[in] type Specified response type.
 ***************************************************************************/
GException::rsp_invalid_type::rsp_invalid_type(string origin,
                                               string type)
{
    m_origin  = origin;
    m_message = "Invalid response type '"+type+"' specified";
}


/***********************************************************************//**
 * @brief General Healpix error
 *
 * @param[in] origin Method that throws the error.
 * @param[in] message Error message.
 ***************************************************************************/
GException::healpix::healpix(string origin, string message)
{
    m_origin  = origin;
    m_message = message;
}


/***********************************************************************//**
 * @brief Healpix resolution invalid
 *
 * @param[in] origin Method that throws the error.
 * @param[in] nside Specified nside parameter.
 ***************************************************************************/
GException::healpix_bad_nside::healpix_bad_nside(string origin, int nside)
{
    m_origin  = origin;
    m_message = "Invalid nside "+str(nside)+" parameter (should be one of "
                "1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192).";
}


/***********************************************************************//**
 * @brief Healpix scheme invalid
 *
 * @param[in] origin Method that throws the error.
 * @param[in] scheme Specified scheme.
 ***************************************************************************/
GException::healpix_bad_scheme::healpix_bad_scheme(string origin, string scheme)
{
    m_origin  = origin;
    m_message = "Invalid scheme '"+scheme+"' (should be either 'RING' or 'NESTED').";
}


/***********************************************************************//**
 * @brief Healpix coordinate system invalid
 *
 * @param[in] origin Method that throws the error.
 * @param[in] coordsys Specified coordinate system.
 ***************************************************************************/
GException::healpix_bad_coords::healpix_bad_coords(string origin, string coordsys)
{
    m_origin  = origin;
    m_message = "Invalid coordinate system '"+coordsys +
                "' (should be either 'EQU' or 'GAL').";
}
