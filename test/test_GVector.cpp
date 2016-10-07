/***************************************************************************
 *                   test_GVector.cpp - Test vector class                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2006-2016 by Juergen Knoedlseder                         *
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
 * @file test_GVector.cpp
 * @brief Testing of vector class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#include <cmath>
#include "test_GVector.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Globals ____________________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/***********************************************************************//**
 * @brief Set parameters and tests
 ***************************************************************************/
void TestGVector::set(void){

    // Test name
    name("GVector");

    // Set parameters
    m_num = 5;

    // Define vectors
    define_vectors();

    // Append tests
    append(static_cast<pfunction>(&TestGVector::allocation), "Vector allocation");
    append(static_cast<pfunction>(&TestGVector::assign), "Assign values");
    append(static_cast<pfunction>(&TestGVector::arithmetics), "Assignment and arithmetics");
    append(static_cast<pfunction>(&TestGVector::comparison), "Comparison");

    return;
}


/***********************************************************************//**
 * @brief Clone test suite
 *
 * @return Pointer to deep copy of test suite.
 ***************************************************************************/
TestGVector* TestGVector::clone(void) const
{
    // Clone test suite
    return new TestGVector(*this);
}


/***********************************************************************//**
 * @brief Define vectors
 ***************************************************************************/
void TestGVector::define_vectors(void){

    // Allocate and initialise vectors
    m_test    = GVector(m_num);
    m_result  = GVector(m_num);
    m_smaller = GVector(m_num+1);
    m_bigger  = GVector(m_num+1);
    for (int i = 0; i < m_num; ++i) {
        m_test[i] = (i+1) * 1.1;
    }
    for (int i = 0; i < m_num-1; ++i) {
        m_smaller[i] = (i+1) * 1.1;
    }
    for (int i = 0; i < m_num+1; ++i) {
        m_bigger[i] = (i+1) * 1.1;
    }
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Vector allocation
 ***************************************************************************/
void TestGVector::allocation(void){

    // Test void constructor
    test_try("Void constructor");
    try {
        GVector vector;
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test zero vector allocation
    test_try("Number constructor");
    try {
        GVector vector(0);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test non-zero allocation
    test_try("Number constructor");
    try {
        GVector vector(10);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test value constructors
    test_try("Value constructors");
    try {
        // One element
        GVector vector1(3.0);
        test_value(vector1[0], 3.0);
        test_assert(vector1.size() == 1, "Expected vector size 1.");

        // Two elements
        GVector vector2(2.0, 5.0);
        test_value(vector2[0], 2.0);
        test_value(vector2[1], 5.0);
        test_assert(vector2.size() == 2, "Expected vector size 2.");

        // Three elements
        GVector vector3(7.0, 8.0, 9.0);
        test_value(vector3[0], 7.0);
        test_value(vector3[1], 8.0);
        test_value(vector3[2], 9.0);
        test_assert(vector3.size() == 3, "Expected vector size 3.");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test copy constructor
    test_try("Copy constructor");
    try {
        GVector vector = m_test;
        test_assert(vector == m_test, vector.print() + " instead of " + m_test.print());
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Vector assignment
 ***************************************************************************/
void TestGVector::assign(void){

    // Test value assignment
    GVector test3(3);
    double  ref = std::acos(-1.0);
    test3[1]    = ref;
    test_value(test3[0], 0.0);
    test_value(test3[1], ref);
    test_value(test3[2], 0.0);
    test_assert(test3.size() == 3, "Expected vector size 3.");

    #if defined(G_RANGE_CHECK)
    test_try("Test out of range access");
    try {
        test3[3];

        test_try_failure();
    }
    catch (GException::out_of_range &e) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }
    #endif

    // Test vector assignment
    m_result = m_test;
    test_assert(m_result == m_test, "m_result == m_test");
    m_result = m_bigger;    
    test_assert(m_result==m_bigger,"m_result == m_bigger");
    test_assert(m_result.size()==m_bigger.size(),"m_result.size() == m_bigger.size()");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Assignment and arithmetics
 ***************************************************************************/
void TestGVector::arithmetics(void){

    m_result  = m_test;
    m_result += m_test;

    // GVector += GVector
    test_assert(m_result[0]==m_test[0]*2&&
                m_result[1]==m_test[1]*2&&
                m_result[2]==m_test[2]*2&&
                m_result[3]==m_test[3]*2&&
                m_result[4]==m_test[4]*2,"GVector += GVector");

    // GVector += 2.0
    m_result  = m_test;
    m_result += 2.0;
    test_value(m_result[0],m_test[0]+2,1e-6,"GVector += 2.0 (1)");
    test_value(m_result[0],m_test[0]+2,1e-6,"GVector += 2.0 (2)");
    test_value(m_result[0],m_test[0]+2,1e-6,"GVector += 2.0 (3)");
    test_value(m_result[0],m_test[0]+2,1e-6,"GVector += 2.0 (4)");

    // GVector -= GVector
    m_result  = m_test;
    m_result -= m_test;
    test_assert(m_result[0]==0&&
                m_result[1]==0&&
                m_result[2]==0&&
                m_result[3]==0&&
                m_result[4]==0,"GVector -= GVector");

    // GVector -= 2.0
    m_result  = m_test;
    m_result -= 2.0;
    test_assert(m_result[0]==m_test[0]-2&&
                m_result[1]==m_test[1]-2&&
                m_result[2]==m_test[2]-2&&
                m_result[3]==m_test[3]-2&&
                m_result[4]==m_test[4]-2,"GVector -= 2.0");

    // GVector *= 2.0
    m_result  = m_test;
    m_result *= 2.0;
    test_assert(m_result[0]==m_test[0]*2&&
                m_result[1]==m_test[1]*2&&
                m_result[2]==m_test[2]*2&&
                m_result[3]==m_test[3]*2&&
                m_result[4]==m_test[4]*2,"GVector *= 2.0");

    // GVector /= 2.0
    m_result  = m_test;
    m_result /= 2.0;
    test_assert(m_result[0]==m_test[0]/2&&
                m_result[1]==m_test[1]/2&&
                m_result[2]==m_test[2]/2&&
                m_result[3]==m_test[3]/2&&
                m_result[4]==m_test[4]/2,"GVector /= 2.0");

    // GVector = -GVector
    m_result = -m_test;
    test_assert(m_result[0]== -m_test[0]&&
                m_result[1]==-m_test[1]&&
                m_result[2]==-m_test[2]&&
                m_result[3]==-m_test[3]&&
                m_result[4]==-m_test[4],"GVector = -GVector");

    // Divide by zero
    m_result  = m_test;
    m_result /= 0.0;
    test_assert(m_result.print() == "(inf, inf, inf, inf, inf)" ||
                m_result.print() == "(Inf, Inf, Inf, Inf, Inf)" ||
                m_result.print() == "(Infinity, Infinity, Infinity, Infinity, Infinity)" ,
                "Divide by zero", m_result.print());

    // GVector + GVector
    m_result = m_test + m_test;
    test_assert(m_result[0]==m_test[0]*2&&
                m_result[1]==m_test[1]*2&&
                m_result[2]==m_test[2]*2&&
                m_result[3]==m_test[3]*2&&
                m_result[4]==m_test[4]*2,"GVector + GVector");

    // GVector + 2.0
    m_result = m_test + 2.0;
    test_value(m_result[0],m_test[0]+2,1e-6,"GVector + 2.0 (1)");
    test_value(m_result[0],m_test[0]+2,1e-6,"GVector + 2.0 (2)");
    test_value(m_result[0],m_test[0]+2,1e-6,"GVector + 2.0 (3)");
    test_value(m_result[0],m_test[0]+2,1e-6,"GVector + 2.0 (4)");

    // GVector + 2
    m_result = m_test + 2;
    test_value(m_result[0],m_test[0]+2,1e-6,"GVector + 2 (1)");
    test_value(m_result[0],m_test[0]+2,1e-6,"GVector + 2 (2)");
    test_value(m_result[0],m_test[0]+2,1e-6,"GVector + 2 (3)");
    test_value(m_result[0],m_test[0]+2,1e-6,"GVector + 2 (4)");

    //
    m_result = 2.0 + m_test;
    test_value(m_result[0],m_test[0]+2,1e-6,"2.0 + GVector (1)");
    test_value(m_result[0],m_test[0]+2,1e-6,"2.0 + GVector (2)");
    test_value(m_result[0],m_test[0]+2,1e-6,"2.0 + GVector (3)");
    test_value(m_result[0],m_test[0]+2,1e-6,"2.0 + GVector (4)");

    // 2 + GVector
    m_result = 2 + m_test;
    test_value(m_result[0],m_test[0]+2,1e-6,"2 + GVector (1)");
    test_value(m_result[0],m_test[0]+2,1e-6,"2 + GVector (2)");
    test_value(m_result[0],m_test[0]+2,1e-6,"2 + GVector (3)");
    test_value(m_result[0],m_test[0]+2,1e-6,"2 + GVector (4)");

    // GVector - GVector
    m_result = m_test - m_test;
    test_assert(m_result[0]==0&&
                m_result[1]==0&&
                m_result[2]==0&&
                m_result[3]==0&&
                m_result[4]==0,"GVector - GVector");

    // GVector - 2.0
    m_result = m_test - 2.0;
    test_assert(m_result[0]==m_test[0]-2&&
                m_result[1]==m_test[1]-2&&
                m_result[2]==m_test[2]-2&&
                m_result[3]==m_test[3]-2&&
                m_result[4]==m_test[4]-2,"GVector - 2.0");

    // 2.0 - GVector
    m_result = 2.0 - m_test;
    test_assert(m_result[0]==2-m_test[0]&&
                m_result[1]==2-m_test[1]&&
                m_result[2]==2-m_test[2]&&
                m_result[3]==2-m_test[3]&&
                m_result[4]==2-m_test[4],"2.0 - GVector");

    // Scalar (or dot) product GVector * GVector
    test_value( m_test[0]*m_test[0]+
                m_test[1]*m_test[1]+
                m_test[2]*m_test[2]+
                m_test[3]*m_test[3]+
                m_test[4]*m_test[4],
                m_test * m_test,
                1e-6,"Scalar (or dot) product GVector * GVector");

    // GVector * 2.0
    m_result = m_test * 2.0;
    test_assert(m_result[0]==m_test[0]*2&&
                m_result[1]==m_test[1]*2&&
                m_result[2]==m_test[2]*2&&
                m_result[3]==m_test[3]*2&&
                m_result[4]==m_test[4]*2,"GVector * 2.0");

    // 2.0 * GVector
    m_result = 2.0 * m_test;
    test_assert(m_result[0]==m_test[0]*2&&
                m_result[1]==m_test[1]*2&&
                m_result[2]==m_test[2]*2&&
                m_result[3]==m_test[3]*2&&
                m_result[4]==m_test[4]*2,"2.0 * GVector");

    // |GVector| (vector norm))
    test_value(sqrt(    m_test[0]*m_test[0]+
                        m_test[1]*m_test[1]+
                        m_test[2]*m_test[2]+
                        m_test[3]*m_test[3]+
                        m_test[4]*m_test[4]),
                        norm(m_test),
                        1e-6,
                        "|GVector| (vector norm)");

    // min(GVector)
    test_value(min(m_test), 1.1, 1e-6, "min(GVector)");

    // max(GVector)
    test_value(max(m_test), 5.5, 1e-6, "max(GVector)");

    // sum(GVector)
    test_value(m_test[0]+m_test[1]+m_test[2]+m_test[3]+m_test[4],
               sum(m_test), 1e-6, "sum(GVector)");

    // acos(GVector/10.0)
    GVector test = acos(m_test/10.0);
    for (int i = 0; i < 5; ++i) {
        test_value(test[i], std::acos(m_test[i]/10.0));
    }

    // acosh(GVector)
    test = acosh(m_test);
    for (int i = 0; i < 5; ++i) {
        test_value(test[i], acosh(m_test[i]));
    }

    // asin(GVector/10.0)
    test = asin(m_test/10.0);
    for (int i = 0; i < 5; ++i) {
        test_value(test[i], std::asin(m_test[i]/10.0));
    }

    // asinh(GVector/10.0)
    test = asinh(m_test/10.0);
    for (int i = 0; i < 5; ++i) {
        test_value(test[i], asinh(m_test[i]/10.0));
    }

    // atan(GVector/10.0)
    test = atan(m_test/10.0);
    for (int i = 0; i < 5; ++i) {
        test_value(test[i], std::atan(m_test[i]/10.0));
    }

    // atanh(GVector/10.0)
    test = atanh(m_test/10.0);
    for (int i = 0; i < 5; ++i) {
        test_value(test[i], atanh(m_test[i]/10.0));
    }

    // cos(GVector)
    test = cos(m_test);
    for (int i = 0; i < 5; ++i) {
        test_value(test[i], std::cos(m_test[i]));
    }

    // cosh(GVector)
    test = cosh(m_test);
    for (int i = 0; i < 5; ++i) {
        test_value(test[i], std::cosh(m_test[i]));
    }

    // exp(GVector)
    test = exp(m_test);
    for (int i = 0; i < 5; ++i) {
        test_value(test[i], std::exp(m_test[i]));
    }

    // abs(cos(m_test))
    test = abs(cos(m_test));
    for (int i = 0; i < 5; ++i) {
        test_value(test[i], std::abs(std::cos(m_test[i])));
    }

    // log(GVector)
    test = log(m_test);
    for (int i = 0; i < 5; ++i) {
        test_value(test[i], std::log(m_test[i]));
    }

    // log10(GVector)
    test = log10(m_test);
    for (int i = 0; i < 5; ++i) {
        test_value(test[i], std::log10(m_test[i]));
    }

    // sin(GVector)
    test = sin(m_test);
    for (int i = 0; i < 5; ++i) {
        test_value(test[i], std::sin(m_test[i]));
    }

    // sinh(GVector)
    test = sinh(m_test);
    for (int i = 0; i < 5; ++i) {
        test_value(test[i], std::sinh(m_test[i]));
    }

    // sqrt(GVector)
    test = sqrt(m_test);
    for (int i = 0; i < 5; ++i) {
        test_value(test[i], std::sqrt(m_test[i]));
    }

    // tan(GVector)
    test = tan(m_test);
    for (int i = 0; i < 5; ++i) {
        test_value(test[i], std::tan(m_test[i]));
    }

    // tanh(GVector)
    test = tanh(m_test);
    for (int i = 0; i < 5; ++i) {
        test_value(test[i], std::tanh(m_test[i]));
    }

    // Incompatible size GVector + GVector
    test_try("Incompatible size GVector + GVector:");
    try {
        m_result = m_test + m_bigger;
        test_try_failure();
    }
    catch (GException::vector_mismatch &e) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // cross(a,b) (using 5-dim vectors)
    test_try("cross(a,b) (using 5-dim vectors)");
    try {
        cross(m_test,m_test);
        test_try_failure();
    }
    catch (GException::vector_bad_cross_dim &e) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    //cross(a,b) (using vectors with different dimension)
    test_try("cross(a,b) (using vectors with different dimension)");
    try {
        cross(m_test,m_bigger);
        test_try_failure();
    }
    catch (GException::vector_mismatch &e) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    //cross(a,b) (using 3-dim vectors)
    test_try("cross(a,b) (using 3-dim vectors)");
    try {
        GVector test_cross_a(3);
        GVector test_cross_b(3);
        test_cross_a[0] = 1.0;
        test_cross_b[1] = 1.0;

        //Test if cross == (0,0,1)
        test_assert(cross(test_cross_a, test_cross_b)[0] == 0 &&
                    cross(test_cross_a, test_cross_b)[1] == 0 &&
                    cross(test_cross_a, test_cross_b)[2] == 1,
                    "Check cross(a,b) value");

        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }
}


/***********************************************************************//**
 * @brief Comparison
 ***************************************************************************/
void TestGVector::comparison(void){

    // GVector == GVector
    test_assert((m_test == m_test),"GVector == GVector");

    // GVector == GVector(0)
    GVector test7(m_num);
    test_assert(!(m_test == test7),"GVector == GVector(0)");

    //GVector == GVector (m_bigger)
    test_assert(!(m_test == m_bigger),"GVector == GVector (m_bigger)");

    // GVector != GVector
    test_assert(!(m_test != m_test),"GVector != GVector");

    // GVector != GVector(0)
    test_assert((m_test != test7),"GVector != GVector(0)");

    // GVector != GVector (m_bigger)
    test_assert((m_test != m_bigger),"GVector != GVector (m_bigger)");
}


/***********************************************************************//**
 * @brief Main test entry point
 ***************************************************************************/
int main(void)
{
    // Allocate test suite container
    GTestSuites testsuites("GVector");

    // Initially assume that we pass all tests
    bool was_successful=true;

    // Create a test suite
    TestGVector test;

    // Append to the container
    testsuites.append(test);

    // Run
    was_successful=testsuites.run();

    // Save xml report
    testsuites.save("reports/GVector.xml");

    // Return
    return was_successful ? 0:1;
}
