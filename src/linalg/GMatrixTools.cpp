/***************************************************************************
 *                     GMatrixTools.cpp  -  Matrix tools                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009 by Jurgen Knodlseder                                *
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
#include "GMatrixTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_SYM_1                       "GSymMatrix sym_matrix(const GMatrix&)"
#define G_SYM_2                 "GSymMatrix sym_matrix(const GMatrixSparse&)"


/***********************************************************************//**
 * @brief GSymMatrix to GMatrix storage class convertor
 *
 * @param[in] m GSymMatrix to be converted
 ***************************************************************************/
GMatrix matrix(const GSymMatrix& m)
{
    // Allocate matrix
    GMatrix result = GMatrix(m.rows(), m.cols());

    // Fill matrix
    for (int col = 0; col < m.cols(); ++col) {
        for (int row = 0; row < m.rows(); ++row)
            result(row, col) = m(row, col);
    }

    // Return
    return result;
}


/***********************************************************************//**
 * @brief GSparseMatrix to GMatrix storage class convertor
 *
 * @param[in] m GSparseMatrix to be converted
 ***************************************************************************/
GMatrix matrix(const GSparseMatrix& m)
{
    // Allocate matrix
    GMatrix result = GMatrix(m.rows(), m.cols());

    // Fill matrix
    for (int col = 0; col < m.cols(); ++col) {
        for (int row = 0; row < m.rows(); ++row)
            result(row, col) = m(row, col);
    }

    // Return
    return result;
}


/***********************************************************************//**
 * @brief GMatrix to GSymMatrix storage class convertor
 *
 * @param[in] m GMatrix to be converted
 ***************************************************************************/
GSymMatrix sym_matrix(const GMatrix& m)
{
    // Allocate matrix
    GSymMatrix result = GSymMatrix(m.rows(), m.cols());

    // Fill matrix
    for (int col = 0; col < m.cols(); ++col) {
        for (int row = col; row < m.rows(); ++row) {
            double value_ll = m(row,col);
            double value_ur = m(col,row);
            if (value_ll != value_ur)
                throw GException::matrix_not_symmetric(G_SYM_1, m.rows(), m.cols());
            result(row, col) = m(row, col);
        }
    }

    // Return
    return result;
}


/***********************************************************************//**
 * @brief GSparseMatrix to GSymMatrix storage class convertor
 *
 * @param[in] m GSparseMatrix to be converted
 ***************************************************************************/
GSymMatrix sym_matrix(const GSparseMatrix& m)
{
    // Allocate matrix
    GSymMatrix result = GSymMatrix(m.rows(), m.cols());

    // Fill matrix
    for (int col = 0; col < m.cols(); ++col) {
        for (int row = col; row < m.rows(); ++row) {
            double value_ll = m(row,col);
            double value_ur = m(col,row);
            if (value_ll != value_ur)
                throw GException::matrix_not_symmetric(G_SYM_2, m.rows(), m.cols());
            result(row, col) = m(row, col);
        }
    }

    // Return
    return result;
}


/***********************************************************************//**
 * @brief GMatrix to GSparseMatrix storage class convertor
 *
 * @param[in] m GMatrix to be converted
 ***************************************************************************/
GSparseMatrix sparse_matrix(const GMatrix& m)
{
    // Allocate matrix
    GSparseMatrix result = GSparseMatrix(m.rows(), m.cols());

    // Fill matrix
    for (int col = 0; col < m.cols(); ++col) {
        for (int row = 0; row < m.rows(); ++row)
            result(row, col) = m(row, col);
    }

    // Return
    return result;
}


/***********************************************************************//**
 * @brief GSymMatrix to GSparseMatrix storage class convertor
 *
 * @param[in] m GSymMatrix to be converted
 ***************************************************************************/
GSparseMatrix sparse_matrix(const GSymMatrix& m)
{
    // Allocate matrix
    GSparseMatrix result = GSparseMatrix(m.rows(), m.cols());

    // Fill matrix
    for (int col = 0; col < m.cols(); ++col) {
        for (int row = 0; row < m.rows(); ++row)
            result(row, col) = m(row, col);
    }

    // Return
    return result;
}

