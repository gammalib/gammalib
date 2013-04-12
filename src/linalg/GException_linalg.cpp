/***************************************************************************
 *            GException_linalg.cpp  -  linalg exception handlers          *
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
 * @file GException_linalg.cpp
 * @brief Exception handlers for linear algebra module
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#include "GException.hpp"
#include "GTools.hpp"


/***********************************************************************//**
 * @brief Empty object exception
 *
 * @param[in] origin Method throwing the exception.
 ***************************************************************************/
GException::empty::empty(std::string origin)
{
    // Set origin
    m_origin  = origin;

    // Set message
    m_message = "Zero-size allocation.";

    // Return
    return;
}


/***********************************************************************//**
 * @brief Index is out of range
 *
 * @param[in] origin Method throwing the exception.
 * @param[in] inx Index.
 * @param[in] min Minimum of valid range.
 * @param[in] max Maximum of valid range.
 ***************************************************************************/
GException::out_of_range::out_of_range(std::string origin,
                                       int         inx,
                                       int         min,
                                       int         max)
{
    // Set origin
    m_origin  = origin;

    // Set message
    m_message = "Index " + gammalib::str(inx) + " out of range [" + gammalib::str(min) +
                "," + gammalib::str(max) + "].";

    // Return
    return;
}


/***********************************************************************//**
 * @brief Value is out of range
 *
 * @param[in] origin Method throwing the exception.
 * @param[in] value Value.
 * @param[in] min Minimum of valid range.
 * @param[in] max Maximum of valid range.
 ***************************************************************************/
GException::out_of_range::out_of_range(std::string origin,
                                       double      value,
                                       double      min,
                                       double      max)
{
    // Set origin
    m_origin  = origin;

    // Set message
    m_message = "Value " + gammalib::str(value) + " out of range [" + gammalib::str(min) +
                "," + gammalib::str(max) + "].";

    // Return
    return;
}


/***********************************************************************//**
 * @brief Vector index is out of range
 *
 * @param[in] origin Method throwing the exception.
 * @param[in] inx Index.
 * @param[in] elements Number of vector elements.
 ***************************************************************************/
GException::out_of_range::out_of_range(std::string origin,
                                       int         inx,
                                       int         elements)
{
    // Set origin
    m_origin = origin;

    // Set message
    if (elements > 0) {
        m_message = "Vector index " + gammalib::str(inx) + " out of range [0," +
                    gammalib::str(elements-1) + "].";
    }
    else {
        m_message = "Empty vector cannot be indexed.";
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Matrix index is out of range
 *
 * @param[in] origin Method throwing the exception.
 * @param[in] row Row index.
 * @param[in] col Column index.
 * @param[in] rows Number of rows in matrix.
 * @param[in] cols Number of columns in matrix.
 ***************************************************************************/
GException::out_of_range::out_of_range(std::string origin,
                                       int         row,
                                       int         col,
                                       int         rows,
                                       int         cols)
{
    // Set origin
    m_origin = origin;

    // Set message
    if (rows > 0 && cols > 0) {
        m_message = "Matrix element (" + gammalib::str(row) + "," + gammalib::str(col) +
                    ") out of range ([0," + gammalib::str(rows-1) + "], [0," +
                    gammalib::str(cols-1) + "])";
    }
    else {
        m_message = "Empty matrix cannot be indexed.";
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Vector dimensions differ
 *
 * @param[in] origin Method throwing the exception.
 * @param[in] size1 Size of first vector.
 * @param[in] size2 Size of second vector.
 ***************************************************************************/
GException::vector_mismatch::vector_mismatch(std::string origin, int size1,
                                             int size2)
{
    // Set origin
    m_origin  = origin;

    // Set message
    m_message = "Vector dimensions differ (" + gammalib::str(size1) + " <-> " +
                gammalib::str(size2) + ").";

    // Return
    return;
}


/***********************************************************************//**
 * @brief Invalid vector dimension for cross product
 *
 * @param[in] origin Method throwing the exception.
 * @param[in] elements Size of vector.
 ***************************************************************************/
GException::vector_bad_cross_dim::vector_bad_cross_dim(std::string origin,
                                                       int elements)
{
    // Set origin
    m_origin  = origin;

    // Set message
    m_message = "Vector cross product only defined for 3 dimensions but"
                " vector size is " + gammalib::str(elements) +".";

    // Return
    return;
}


/***********************************************************************//**
 * @brief Mismatch between matrix and vector
 *
 * @param[in] origin Method throwing the exception.
 * @param[in] num Vector length.
 * @param[in] rows Number of matrix rows.
 * @param[in] cols Number of matrix columns.
 ***************************************************************************/
GException::matrix_vector_mismatch::matrix_vector_mismatch(std::string origin,
                                                           int         num, 
                                                           int         rows,
                                                           int         cols)
{
    // Set origin
    m_origin  = origin;

    // Set message
    m_message = "Vector length " + gammalib::str(num) +
                " is incompatible with matrix size (" +
                gammalib::str(rows) + "," + gammalib::str(cols) + ").";

    // Return
    return;
}


/***********************************************************************//**
 * @brief Mismatch between matrices
 *
 * @param[in] origin Method throwing the exception.
 * @param[in] rows1 Number of rows in first matrix.
 * @param[in] cols1 Number of columns in first matrix.
 * @param[in] rows2 Number of rows in second matrix.
 * @param[in] cols2 Number of columns in second matrix.
 ***************************************************************************/
GException::matrix_mismatch::matrix_mismatch(std::string origin, int rows1,
                                             int cols1, int rows2, int cols2)
{
    // Set origin
    m_origin  = origin;

    // Set message
    m_message  = "Matrices are incompatible for operation.";
    m_message += "Matrix size ";
    m_message += "(" + gammalib::str(rows1) + "," + gammalib::str(cols1) + ")";
    m_message += " differs from matrix size ";
    m_message += "(" + gammalib::str(rows2) + "," + gammalib::str(cols2) + ").";

    // Return
    return;
}


/***********************************************************************//**
 * @brief Matrix not square
 *
 * @param[in] origin Method throwing the exception.
 * @param[in] rows Number of matrix rows.
 * @param[in] cols Number of matrix columns.
 ***************************************************************************/
GException::matrix_not_square::matrix_not_square(std::string origin,
                                                 int         rows,
                                                 int         cols)
{
    // Set origin
    m_origin  = origin;

    // Set message
    m_message  = "Matrix of size (" + gammalib::str(rows) + "," + gammalib::str(cols) + ")";
    m_message += " is not square.";

    // Return
    return;
}


/***********************************************************************//**
 * @brief Matrix is not positive definite
 *
 * @param[in] origin Method throwing the exception.
 * @param[in] row Matrix row.
 * @param[in] sum Row sum.
 ***************************************************************************/
GException::matrix_not_pos_definite::matrix_not_pos_definite(std::string origin,
                                                             int row, double sum)
{
    // Set origin
    m_origin  = origin;

    // Set message
    m_message = "Matrix is not positive definite (sum " + gammalib::str(sum) +
                " occured in row/column " + gammalib::str(row) + ").";

    // Return
    return;
}


/***********************************************************************//**
 * @brief Matrix is not symmetric
 *
 * @param[in] origin Method throwing the exception.
 * @param[in] rows Number of rows.
 * @param[in] cols Number of columns.
 ***************************************************************************/
GException::matrix_not_symmetric::matrix_not_symmetric(std::string origin,
                                                       int rows, int cols)
{
    // Set origin
    m_origin  = origin;

    // Set message
    m_message = "Matrix (" + gammalib::str(rows) + "," + gammalib::str(cols) + ") is not symmetric.";

    // Return
    return;
}


/***********************************************************************//**
 * @brief Matrix has not been factorised
 *
 * @param[in] origin Method throwing the exception.
 * @param[in] type Factorisation type.
 ***************************************************************************/
GException::matrix_not_factorised::matrix_not_factorised(std::string origin,
                                                         std::string type)
{
    // Set origin
    m_origin  = origin;

    // Set message
    m_message = "Matrix has not been factorised using " + type;

    // Return
    return;
}


/***********************************************************************//**
 * @brief All matrix elements are zero
 *
 * @param[in] origin Method throwing the exception.
 ***************************************************************************/
GException::matrix_zero::matrix_zero(std::string origin)
{
    // Set origin
    m_origin = origin;

    // Set message
    m_message = "All matrix elements are zero.";

    // Return
    return;
}


/***********************************************************************//**
 * @brief Invalid ordering scheme requested
 *
 * @param[in] origin Method throwing the exception.
 * @param[in] order Ordering scheme.
 * @param[in] min_order Minimum ordering scheme.
 * @param[in] max_order Maximum ordering scheme.
 ***************************************************************************/
GException::invalid_order::invalid_order(std::string origin, int order,
                                         int min_order, int max_order)
{
    // Set origin
    m_origin  = origin;

    // Set message
    m_message = "Invalid ordering type " + gammalib::str(order) + 
                "requested; must be comprised in [" + gammalib::str(min_order) +
                "," + gammalib::str(max_order) + "]";

    // Return
    return;
}
