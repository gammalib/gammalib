/***************************************************************************
 *              test_GMatrix.hpp - Test generic matrix class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2021 by Juergen Knoedlseder                         *
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
 * @file test_GMatrix.hpp
 * @brief Definition of unit tests for generic matrices
 * @author Juergen Knoedlseder
 */

#ifndef TEST_GMATRIX_HPP
#define TEST_GMATRIX_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"


/***********************************************************************//**
 * @class TestGMatrix
 *
 * @brief Test suite for generic matrix class testing
 *
 * This test suite performs a unit test for the generic matrix class. Tests
 * that match logically are gathered in a number of class methods, and the
 * class methods are executed when running the test.
 ***************************************************************************/
class TestGMatrix : public GTestSuite {

public:
    // Constructors and destructors
    TestGMatrix(void) : GTestSuite() {}
    virtual ~TestGMatrix(void) {}

    // Methods
    virtual void         set(void);
    virtual TestGMatrix* clone(void) const;
    virtual std::string  classname(void) const { return "TestGMatrix"; }
    void                 empty(void);
    void                 alloc_matrix(void);
    void                 assign_values(void);
    void                 copy_matrix(void);
    void                 matrix_operations(void);
    void                 matrix_arithmetics(void);
    void                 matrix_functions(void);
    void                 matrix_compare(void);
    void                 matrix_print(void);

private:
    // Private methods
    GMatrix set_matrix(void) const;
    GMatrix set_matrix_zero(void) const;
    GVector set_vector(void) const;
    bool    check_matrix(const GMatrix& matrix,
                         const double&     scale = 1.0,
                         const double&     offset = 0.0) const;
    bool    check_matrix_trans(const GMatrix& matrix,
                               const double&  scale,
                               const double&  offset) const;
    bool    check_matrix_lt(const GMatrix& matrix, const GMatrix& ref) const;
    bool    check_matrix_ut(const GMatrix& matrix, const GMatrix& ref) const;
    
    // Private members;
    GMatrix m_test;    //!< Test matrix
    GMatrix m_bigger;  //!< Bigger test matrix (for collisions)
    GVector v_test;    //!< Test vector
};

#endif /* TEST_GMATRIX_HPP */
