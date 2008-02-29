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


/***************************************************************************
 *                               FITS error                                *
 ***************************************************************************/
GException::fits_error::fits_error(string origin, int status)
{
    m_origin  = origin;
    ostringstream s_error;
    char err_text[31];
    __ffgerr(status, err_text);
    s_error << err_text << " (status: " << status << ")";
    m_message = s_error.str();
}


/***************************************************************************
 *                          FITS file opening error                        *
 ***************************************************************************/
GException::fits_open_error::fits_open_error(string origin, string filename, int status)
{
    m_origin  = origin;
    m_message = "Unable to open FITS file '" + filename + "'";
}


/***************************************************************************
 *                          FITS file already opened                       *
 ***************************************************************************/
GException::fits_already_opened::fits_already_opened(string origin, string filename)
{
    m_origin  = origin;
    m_message = "GFits object has already the FITS file '" + filename + "' opened";
}


/***************************************************************************
 *                        FITS keyword not found error                     *
 ***************************************************************************/
GException::fits_key_not_found::fits_key_not_found(string origin, string keyname, int status)
{
    m_origin  = origin;
    m_message = "Keyword '" + keyname + "' not found in header";
}


/***************************************************************************
 *                         FITS column not found error                     *
 ***************************************************************************/
GException::fits_column_not_found::fits_column_not_found(string origin, string colname, int status)
{
    m_origin  = origin;
    m_message = "Column '" + colname + "' not found in table";
}


/***************************************************************************
 *                           FITS HDU not found error                      *
 ***************************************************************************/
GException::fits_hdu_not_found::fits_hdu_not_found(string origin, string extname, int status)
{
    m_origin  = origin;
    m_message = "HDU '" + extname + "' not found in FITS file";
}


/***************************************************************************
 *                          Unknown HDU type found                         *
 ***************************************************************************/
GException::fits_unknown_HDU_type::fits_unknown_HDU_type(string origin, int type)
{
    m_origin  = origin;
    m_message = "HDU type " + str(type) + " is not known";
}


/***************************************************************************
 *                             HDU is not a table                          *
 ***************************************************************************/
GException::fits_HDU_not_a_table::fits_HDU_not_a_table(string origin, int type)
{
    m_origin  = origin;
    m_message = "HDU is not a table (type " + str(type) + ")";
}


/***************************************************************************
 *                             HDU is not a table                          *
 ***************************************************************************/
GException::fits_invalid_type::fits_invalid_type(string origin, string message)
{
    m_origin  = origin;
    m_message = message;
}


/***************************************************************************
 *                             Column type is unknown                      *
 ***************************************************************************/
GException::fits_unknown_coltype::fits_unknown_coltype(string origin, int type)
{
    m_origin  = origin;
    m_message = "Column type " + str(type) + " is unknown";
}
