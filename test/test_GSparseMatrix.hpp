/***************************************************************************
 *            test_GSparseMatrix.hpp  -  Test sparse matrix class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012 by Juergen Knoedlseder                              *
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
 * @file test_GSparseMatrix.hpp
 * @brief Definition of unit tests for sparse matrices
 * @author Juergen Knoedlseder
 */

#ifndef TEST_GSPARSEMATRIX_HPP
#define TEST_GSPARSEMATRIX_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"


/***********************************************************************//**
 * @class TestGSparseMatrix
 *
 * @brief Test suite for sparse matrix class testing
 *
 * This test suite performs a unit test for the sparse matrix class. Tests
 * that match logically are gathered in a number of class methods, and the
 * class methods are executed when running the test.
 ***************************************************************************/
class TestGSparseMatrix : public GTestSuite {

public:
    // Constructors and destructors
    TestGSparseMatrix(void) : GTestSuite() {}
    virtual ~TestGSparseMatrix(void) {}

    // Methods
    virtual void set(void);
    void         alloc_matrix(void);
    void         assign_values(void);
    void         copy_matrix(void);
    void         matrix_operations(void);
    void         matrix_arithmetics(void);
    void         matrix_functions(void);
    void         matrix_compare(void);
    void         matrix_cholesky(void);
    void         matrix_print(void);

private:
    // Private methods
    GSparseMatrix set_matrix(void) const;
    //GSparseMatrix set_matrix_zero(void) const;
    GVector       set_vector(void) const;
    bool          check_matrix(const GSparseMatrix& matrix,
                               const double&     scale = 1.0,
                               const double&     offset = 0.0) const;
    bool          check_matrix_trans(const GSparseMatrix& matrix,
                                     const double&  scale,
                                     const double&  offset) const;
    bool          check_matrix_lt(const GSparseMatrix& matrix,
                                  const GSparseMatrix& ref) const;
    bool          check_matrix_ut(const GSparseMatrix& matrix,
                                  const GSparseMatrix& ref) const;
    
    // Private members;
    GSparseMatrix m_test;    //!< Test matrix
    GSparseMatrix m_bigger;  //!< Bigger test matrix (for collisions)
    GVector       v_test;    //!< Test vector
};

#endif /* TEST_GSPARSEMATRIX_HPP */
