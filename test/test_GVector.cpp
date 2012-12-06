/***************************************************************************
 *                  test_GVector.cpp  -  test vector class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2006-2012 by Jurgen Knodlseder                           *
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
#include <cmath>
#include <iostream>
#include <stdexcept>                          // std::exception
#include "test_GVector.hpp"

void TestGVector::set(void){
    // Test name
    name("GVector");

    // Set parameters
    m_num=5;

    //add tests
    add_test(static_cast<pfunction>(&TestGVector::define_vectors),"Define vectors");
    add_test(static_cast<pfunction>(&TestGVector::test1),"Test 1: Allocate zero vector");
    add_test(static_cast<pfunction>(&TestGVector::test2),"Test 2: Allocate too large vector");
    add_test(static_cast<pfunction>(&TestGVector::test3),"Test 3: Assign values");
    add_test(static_cast<pfunction>(&TestGVector::test4),"Test 4: Define vector using copy constructor");
    add_test(static_cast<pfunction>(&TestGVector::test5),"Test 5: Vector assignment");
    add_test(static_cast<pfunction>(&TestGVector::test6),"Test 6: Assignment and arithmetics");
    add_test(static_cast<pfunction>(&TestGVector::test7),"Test 7: Comparison");

    return;
}

// Define vectors
void TestGVector::define_vectors(void){
    m_test= GVector(m_num);
    m_result=GVector(m_num);
    m_smaller=GVector(m_num+1);
    m_bigger=GVector(m_num+1);
    for (int i = 0; i < m_num; ++i)
        m_test[i] = (i+1) * 1.1;
    for (int i = 0; i < m_num-1; ++i)
        m_smaller[i] = (i+1) * 1.1;
    for (int i = 0; i < m_num+1; ++i)
        m_bigger[i] = (i+1) * 1.1;
}

//Test 1: Allocate zero vector
void TestGVector::test1(void){
    GVector test1(0);
}

//Test 2: Allocate too large vector
void TestGVector::test2(void){
    /*
    GVector test2(1000000000);
    */
}

//Test 3: Assign values
void TestGVector::test3(void){
    GVector test3(3);
    test3[1] = acos(-1.0);

    test_value(test3[1],acos(-1.0),1e-6,"test3 == (0, pi, 0)");
    test_assert(test3.size()==3,"test3.size()==3");

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
}

//Test 4: Vector copy constructor
void TestGVector::test4(void){

    GVector test4 = m_test;
    test_assert(test4==m_test,"test4 == m_test");
}

//Test 5: Vector assignment
void TestGVector::test5(void){
    m_result = m_test;

    test_assert(m_result==m_test,"m_result == m_test");

    m_result = m_bigger;    
    test_assert(m_result==m_bigger,"m_result == m_bigger");
    test_assert(m_result.size()==m_bigger.size(),"m_result.size() == m_bigger.size()");
}

//Test 6: Assignment and arithmetics
void TestGVector::test6(void){
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

    // Devide by zero
    m_result  = m_test;
    m_result /= 0.0;
    test_assert(m_result.print()=="(inf, inf, inf, inf, inf)"||m_result.print()=="(Inf, Inf, Inf, Inf, Inf)","Devide by zero",m_result.print());

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
    test_assert(min(m_test)==1.1,"min(GVector)");

    // max(GVector)
    test_assert(max(m_test)==5.5,"max(GVector)");

    // sum(GVector)
    test_value(m_test[0]+m_test[1]+m_test[2]+m_test[3]+m_test[4],sum(m_test),1e-6,"sum(GVector)");

    // acos(GVector/10.0)
    test_assert(acos(m_test/10.0).print()=="(1.46057, 1.34898, 1.23449, 1.1152, 0.988432)","acos(GVector/10.0)");

    // acosh(GVector)
    test_assert(acosh(m_test).print()=="(0.443568, 1.42542, 1.86328, 2.16158, 2.38953)","acosh(GVector)");

    // asin(GVector/10.0)
    test_assert(asin(m_test/10.0).print()=="(0.110223, 0.221814, 0.336304, 0.455599, 0.582364)","asin(GVector/10.0)");

    // asinh(GVector/10.0)
    test_assert(asinh(m_test/10.0).print()=="(0.109779, 0.218263, 0.324286, 0.426913, 0.52548)","asinh(GVector/10.0)");

    // atan(GVector/10.0)
    test_assert(atan(m_test/10.0).print()=="(0.10956, 0.21655, 0.318748, 0.414507, 0.502843)","atan(GVector/10.0)");

    // atanh(GVector/10.0)
    test_assert(atanh(m_test/10.0).print()=="(0.110447, 0.223656, 0.342828, 0.472231, 0.618381)","atanh(GVector/10.0)");

    // cos(GVector)
    test_assert(cos(m_test).print()=="(0.453596, -0.588501, -0.98748, -0.307333, 0.70867)","cos(GVector)");

    // cosh(GVector)
    test_assert(cosh(m_test).print()=="(1.66852, 4.56791, 13.5748, 40.7316, 122.348)","cosh(GVector)");

    // exp(GVector)
    test_assert(exp(m_test).print()=="(3.00417, 9.02501, 27.1126, 81.4509, 244.692)","exp(GVector)");

    // abs(cos(m_test))
    test_assert(abs(cos(m_test)).print()=="(0.453596, 0.588501, 0.98748, 0.307333, 0.70867)","abs(cos(m_test))");

    // log(GVector)
    test_assert(log(m_test).print()=="(0.0953102, 0.788457, 1.19392, 1.4816, 1.70475)","log(GVector)");

    // log10(GVector)
    test_assert(log10(m_test).print()=="(0.0413927, 0.342423, 0.518514, 0.643453, 0.740363)","log10(GVector)");

    // sin(GVector)
    test_assert(sin(m_test).print()=="(0.891207, 0.808496, -0.157746, -0.951602, -0.70554)","sin(GVector)");

    // sinh(GVector)
    test_assert(sinh(m_test).print()=="(1.33565, 4.45711, 13.5379, 40.7193, 122.344)","sinh(GVector)");

    // sqrt(GVector)
    test_assert(sqrt(m_test).print()=="(1.04881, 1.48324, 1.81659, 2.09762, 2.34521)","sqrt(GVector)");

    // tan(GVector)
    test_assert(tan(m_test).print()=="(1.96476, -1.37382, 0.159746, 3.09632, -0.995584)","tan(GVector)");

    // tanh(GVector)
    test_assert(tanh(m_test).print()=="(0.800499, 0.975743, 0.997283, 0.999699, 0.999967)","tanh(GVector)");

    //

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
        test_assert(cross(test_cross_a, test_cross_b)[0]==0&&cross(test_cross_a, test_cross_b)[1]==0&&cross(test_cross_a, test_cross_b)[2]==1,"Check cross(a,b) value");

        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }
}

//Test 7: Comparison
void TestGVector::test7(void){

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
int main(void)
{
    GTestSuites testsuites("GVector");

    bool was_successful=true;

    //Create a test suite
    TestGVector test;

    //Append to the container
    testsuites.append(test);

    //Run
    was_successful=testsuites.run();

    //save xml report
    testsuites.save("reports/GVector.xml");

    // Return
    return was_successful ? 0:1;
}
