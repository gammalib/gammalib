/***************************************************************************
 *              test_GVector.hpp  -   test vector class                    *
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

#ifndef TEST_GVECTOR_HPP
#define TEST_GVECTOR_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"
#include <iostream>                           // cout, cerr

class TestGVector : public GTestSuite
{
    public:
        // Constructors and destructors
        TestGVector(void) : GTestSuite(){ return; }
        virtual ~TestGVector(void){ return; }

        // Methods
        virtual void set(void);
        void define_vectors(void);
        void test1(void);
        void test2(void);
        void test3(void);
        void test4(void);
        void test5(void);
        void test6(void);
        void test7(void);


    // Private members
    private:
        int m_num;
        GVector m_test;
        GVector m_result;
        GVector m_smaller;
        GVector m_bigger;
};

#endif /* TEST_GVECTOR_HPP */
