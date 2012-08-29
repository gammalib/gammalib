/***************************************************************************
 *              test_GMatrix.hpp  -   test matrix class                    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012 by Jean-Baptiste Cayrou                             *
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

#ifndef TEST_GMATRIX_HPP
#define TEST_GMATRIX_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"
#include <iostream>                           // cout, cerr

class TestGMatrix : public GTestSuite
{
    public:
        // Constructors and destructors
        TestGMatrix(void) : GTestSuite(){ return; }
        virtual ~TestGMatrix(void){ return; }

        // Methods
        virtual void set(void);

        void set_matrix(void);
        void set_vector(void);
        void set_bigger(void);
        void test_output(void);
        void test_conversion(void);
        void test_extract(void);
        void test1(void);
        void test2(void);
        void test3(void);
        void test4(void);
        void test5(void);
        void test6(void);
        void test7(void);
        void test8(void);
        void test9(void);
        void test10(void);
        void test11(void);

    // Private methods
    private:
        void init_matrix(void);
        void init_vector(void);

        int check_matrix(const GMatrix& m, const double scale, const double add);
        int check_transpose_matrix(const GMatrix& m, const double scale, const double add);
        int check_matrix_lt(const GMatrix& m, const double scale, const double add);
        int check_matrix_ut(const GMatrix& m, const double scale, const double add);
        int check_matrix_vector(const GVector& v);
        int check_matrix_matrix(const GMatrix& m);
        int check_matrix_min(const double min);
        int check_matrix_max(const double max);
        int check_matrix_sum(const double sum);



    
    // Private attributes
        GMatrix  m_test;
        GMatrix  bigger;
        GMatrix  result;
        GVector  v_test;
        double   m_g_matrix[12];//
        double   m_g_vector[4];
        int      m_g_rows;
        int      m_g_cols;

};

#endif /* TEST_GMATRIX_HPP */
