/***************************************************************************
 *             test_GMatrixSparse.hpp - Test sparse matrix class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2013 by Juergen Knoedlseder                         *
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
 * @file test_GMatrixSparse.hpp
 * @brief Definition of unit tests for sparse matrices
 * @author Juergen Knoedlseder
 */

#ifndef TEST_GMATRIXSPARSE_HPP
#define TEST_GMATRIXSPARSE_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"


/***********************************************************************//**
 * @class TestGMatrixSparse
 *
 * @brief Test suite for sparse matrix class testing
 *
 * This test suite performs a unit test for the sparse matrix class. Tests
 * that match logically are gathered in a number of class methods, and the
 * class methods are executed when running the test.
 ***************************************************************************/
class TestGMatrixSparse : public GTestSuite {

public:
    // Constructors and destructors
    TestGMatrixSparse(void) : GTestSuite() {}
    virtual ~TestGMatrixSparse(void) {}

    // Methods
    virtual void               set(void);
    virtual TestGMatrixSparse* clone(void) const;
    void                       alloc_matrix(void);
    void                       assign_values(void);
    void                       copy_matrix(void);
    void                       matrix_operations(void);
    void                       matrix_arithmetics(void);
    void                       matrix_functions(void);
    void                       matrix_compare(void);
    void                       matrix_cholesky(void);
    void                       matrix_print(void);

private:
    // Private methods
    GMatrixSparse set_matrix(void) const;
    //GMatrixSparse set_matrix_zero(void) const;
    GVector       set_vector(void) const;
    bool          check_matrix(const GMatrixSparse& matrix,
                               const double&        scale = 1.0,
                               const double&        offset = 0.0) const;
    bool          check_matrix_trans(const GMatrixSparse& matrix,
                                     const double&        scale,
                                     const double&        offset) const;
    bool          check_matrix_lt(const GMatrixSparse& matrix,
                                  const GMatrixSparse& ref) const;
    bool          check_matrix_ut(const GMatrixSparse& matrix,
                                  const GMatrixSparse& ref) const;
    
    // Private members;
    GMatrixSparse m_test;    //!< Test matrix
    GMatrixSparse m_bigger;  //!< Bigger test matrix (for collisions)
    GVector       v_test;    //!< Test vector
};

#endif /* TEST_GMATRIXSPARSE_HPP */
