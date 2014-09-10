/***************************************************************************
 *                 test_GVector.hpp - Test vector class                    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2014 by Jean-Baptiste Cayrou                        *
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
 * @file test_GVector.hpp
 * @brief Testing of vector class definition
 * @author Juergen Knoedlseder
 */

#ifndef TEST_GVECTOR_HPP
#define TEST_GVECTOR_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"


/***********************************************************************//**
 * @class TestGVector
 *
 * @brief Test suite for GVector class
 ***************************************************************************/
class TestGVector : public GTestSuite
{
public:
    // Constructors and destructors
    TestGVector(void) : GTestSuite(), m_num(0) { }
    virtual ~TestGVector(void) { }

    // Methods
    virtual void         set(void);
    virtual TestGVector* clone(void) const;
    virtual std::string  classname(void) const { return "TestGVector"; }
    void                 define_vectors(void);
    void                 allocation(void);
    void                 assign(void);
    void                 arithmetics(void);
    void                 comparison(void);

// Private members
private:
    int     m_num;
    GVector m_test;
    GVector m_result;
    GVector m_smaller;
    GVector m_bigger;
};

#endif /* TEST_GVECTOR_HPP */
