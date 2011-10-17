/***************************************************************************
 *            GException_linalg.cpp  -  linalg exception handlers          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
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

/* __ Includes ___________________________________________________________ */
#include "GException.hpp"
#include "GTools.hpp"


/***************************************************************************
 *                           Empty object exception                        *
 ***************************************************************************/
GException::empty::empty(std::string origin)
{
    m_origin  = origin;
    m_message = "Zero-size allocation";
}


/***************************************************************************
 *                            Index out of range                           *
 ***************************************************************************/
GException::out_of_range::out_of_range(std::string origin,
                                       int         inx,
                                       int         min,
                                       int         max)
{
    m_origin  = origin;
    m_message = "Index (" + str(inx) + ") out of range [" + str(min) +
                "," + str(max) + "]";
}


/***************************************************************************
 *                            Value out of range                           *
 ***************************************************************************/
GException::out_of_range::out_of_range(std::string origin, double value,
                                       double min, double max)
{
    m_origin  = origin;
    m_message = "Value (" + str(value) + ") out of range [" + str(min) +
                "," + str(max) + "]";
}


/***************************************************************************
 *                          Vector index out of range                      *
 ***************************************************************************/
GException::out_of_range::out_of_range(std::string origin,
                                       int         inx,
                                       int         elements)
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
GException::out_of_range::out_of_range(std::string origin,
                                       int         row,
                                       int         col,
                                       int         rows,
                                       int         cols)
{
    m_origin  = origin;
    m_message = "Matrix element (" + str(row) + "," + str(col) +
                ") out of range ([0," + str(rows) + "], [0," +
                str(cols) + "])";
}


/***************************************************************************
 *                          Vector dimensions differ                       *
 ***************************************************************************/
GException::vector_mismatch::vector_mismatch(std::string origin, int size1,
                                             int size2)
{
    m_origin  = origin;
    m_message = "Vector dimensions differ (" + str(size1) + " <-> " +
                str(size2) + ")";
}


/***************************************************************************
 *                   Invalid vector dimension for cross product            *
 ***************************************************************************/
GException::vector_bad_cross_dim::vector_bad_cross_dim(std::string origin,
                                                       int elements)
{
    m_origin  = origin;
    m_message = "Vector cross product only defined for 3 dimensions but"
                " vector size is " + str(elements);
}


/***************************************************************************
 *                     Mismatch between matrix and vector                  *
 ***************************************************************************/
GException::matrix_vector_mismatch::matrix_vector_mismatch(std::string origin,
                                                           int num, int rows,
                                                           int cols)
{
    m_origin  = origin;
    m_message = "Vector dimension [" + str(num) +
                "] is incompatible with matrix size [" +
                str(rows) + "," + str(cols) + "]";
}


/***************************************************************************
 *                       Matrix mismatch in operation                      *
 ***************************************************************************/
GException::matrix_mismatch::matrix_mismatch(std::string origin, int rows1,
                                             int cols1, int rows2, int cols2)
{
    m_origin  = origin;
    m_message = "Matrix mismatch: M1(" + str(rows1) + "," + str(cols1) +
                ") incompatible with M2(" + str(rows2) + "," + str(cols2) + ")";
}


/***************************************************************************
 *                           Matrix not square                             *
 ***************************************************************************/
GException::matrix_not_square::matrix_not_square(std::string origin,
                                                 int         rows,
                                                 int         cols)
{
    m_origin  = origin;
    m_message = "Matrix is not square [" + str(rows) + "," + str(cols) + "]";
}


/***************************************************************************
 *                     Matrix is not positive definite                     *
 ***************************************************************************/
GException::matrix_not_pos_definite::matrix_not_pos_definite(std::string origin,
                                                             int row, double sum)
{
    m_origin  = origin;
    m_message = "Matrix is not positive definite (sum " + str(sum) +
                " occured in row/column " + str(row) + ")";
}


/***************************************************************************
 *                         Matrix is not symmetric                         *
 ***************************************************************************/
GException::matrix_not_symmetric::matrix_not_symmetric(std::string origin,
                                                       int  cols, int rows)
{
    m_origin  = origin;
    m_message = "Matrix is not symmetric [" + str(rows) + "," + str(cols) + "]";
}


/***************************************************************************
 *                       Matrix has not been factorised                    *
 ***************************************************************************/
GException::matrix_not_factorised::matrix_not_factorised(std::string origin,
                                                         std::string type)
{
    m_origin  = origin;
    m_message = "Matrix has not been factorised using " + type;
}


/***************************************************************************
 *                       All matrix elements are zero                      *
 ***************************************************************************/
GException::matrix_zero::matrix_zero(std::string origin)
{
    m_origin = origin;
    m_message = "All matrix elements are zero";
}


/***************************************************************************
 *                     Invalid ordering scheme requested                   *
 ***************************************************************************/
GException::invalid_order::invalid_order(std::string origin, int order,
                                         int min_order, int max_order)
{
    m_origin  = origin;
    m_message = "Invalid ordering type " + str(order) + 
                "requested; must be comprised in [" + str(min_order) +
                "," + str(max_order) + "]";
}
