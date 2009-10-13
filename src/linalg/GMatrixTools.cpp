/***************************************************************************
 *                     GMatrixTools.cpp  -  Matrix tools                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

/* __ Includes ___________________________________________________________ */
#include "GMatrixTools.hpp"


/***********************************************************************//**
 * @brief GSymMatrix to GMatrix storage class convertor
 *
 * @param[in] m GSymMatrix to be converted
 ***************************************************************************/
GMatrix convert(const GSymMatrix& m)
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
GMatrix convert(const GSparseMatrix& m)
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
