/***************************************************************************
 *               GException_fits.cpp  -  fits exception handlers           *
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

/* __ Includes ___________________________________________________________ */
#include "GException.hpp"
#include "GTools.hpp"
#include "GFitsCfitsio.hpp"


/***********************************************************************//**
 * @brief General FITS error
 *
 * @param[in] origin Method that throws the error.
 * @param[in] status cfitsio status.
 * @param[in] message Optional error message.
 ***************************************************************************/
GException::fits_error::fits_error(std::string origin, 
                                   int         status,
                                   std::string message)
{
    m_origin  = origin;
    std::ostringstream s_error;
    char err_text[31];
    __ffgerr(status, err_text);
    s_error << err_text << " (status=" << status << ") " + message;
    m_message = s_error.str();
}


/***********************************************************************//**
 * @brief FITS error: unable to open FITS file
 *
 * @param[in] origin Method that throws the error.
 * @param[in] filename Name of the file for which opening was attempted.
 * @param[in] status cfitsio status.
 ***************************************************************************/
GException::fits_open_error::fits_open_error(std::string origin,
                                             std::string filename,
                                             int         status)
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
GException::fits_file_exist::fits_file_exist(std::string origin,
                                             std::string filename,
                                             int         status)
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
GException::fits_already_opened::fits_already_opened(std::string origin,
                                                     std::string filename)
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
GException::fits_key_not_found::fits_key_not_found(std::string origin,
                                                   std::string keyname,
                                                   int         status)
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
GException::fits_column_not_found::fits_column_not_found(std::string origin,
                                                         std::string colname,
                                                         int         status)
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
GException::fits_no_header::fits_no_header(std::string origin,
                                           std::string message,
                                           int         status)
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
GException::fits_no_data::fits_no_data(std::string origin,
                                       std::string message,
                                       int         status)
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
GException::fits_hdu_not_found::fits_hdu_not_found(std::string origin,
                                                   std::string extname,
                                                   int         status)
{
    m_origin  = origin;
    m_message = "HDU \"" + extname + "\" not found in FITS file";
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
GException::fits_hdu_not_found::fits_hdu_not_found(std::string origin,
                                                   int         extno,
                                                   int         status)
{
    m_origin  = origin;
    m_message = "HDU number " + str(extno) + " not found in FITS file";
    if (status != 0)
        m_message += " (status=" + str(status) + ")";
}


/***********************************************************************//**
 * @brief FITS error: HDU of unknown type found
 *
 * @param[in] origin Method that throws the error.
 * @param[in] type Specified HDU type.
 ***************************************************************************/
GException::fits_unknown_HDU_type::fits_unknown_HDU_type(std::string origin,
                                                         int         type)
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
GException::fits_HDU_not_a_table::fits_HDU_not_a_table(std::string origin,
                                                       int         type)
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
GException::fits_invalid_type::fits_invalid_type(std::string origin,
                                                 std::string message)
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
GException::fits_unknown_tabtype::fits_unknown_tabtype(std::string origin,
                                                       int         type)
{
    m_origin  = origin;
    m_message = "Table type '" + str(type) + "' is unknown";
}


/***********************************************************************//**
 * @brief FITS error: Column type is unknown
 *
 * @param[in] origin Method that throws the error.
 * @param[in] colname Column name.
 * @param[in] type Specified column type.
 ***************************************************************************/
GException::fits_unknown_coltype::fits_unknown_coltype(std::string origin,
                                                       std::string colname,
                                                       int         type)
{
    m_origin  = origin;
    m_message = "Column \""+colname+"\" has unsupported typecode="+str(type);
}


/***********************************************************************//**
 * @brief FITS error: Bad column length
 *
 * @param[in] origin Method that throws the error.
 * @param[in] length Length of the column.
 * @param[in] rows Number of rows in table.
 ***************************************************************************/
GException::fits_bad_col_length::fits_bad_col_length(std::string origin,
                                                     int         length,
                                                     int         rows)
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
GException::fits_bad_bitpix::fits_bad_bitpix(std::string origin,
                                             int         bitpix)
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
GException::fits_wrong_image_operator::fits_wrong_image_operator(std::string origin,
                                                                 int         naxis,
                                                                 int         nargs)
{
    m_origin  = origin;
    m_message = "Wrong image pixel access operator has been used (dimension=" +
                str(naxis) + " <=> arguments=" + str(nargs) + ")";
}
