/***************************************************************************
 *                test_GNumerics.cpp  -  test numerics modules             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2018 by Jurgen Knodlseder                           *
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
 * @file test_GNumerics.cpp
 * @brief Unit tests implementation for numerics module
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <ostream>
#include <stdexcept>
#include <stdlib.h>
#include "test_GNumerics.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Globals ____________________________________________________________ */


/***********************************************************************//**
 * @brief Set tests
 ***************************************************************************/
void TestGNumerics::set(void)
{
    //set name
    name("GNumerics");

    // Set parameters
    m_sigma = 2.5;

    // Append tests
    append(static_cast<pfunction>(&TestGNumerics::test_math),
           "Test GMath");
    append(static_cast<pfunction>(&TestGNumerics::test_ndarray),
           "Test GNdarray");
    append(static_cast<pfunction>(&TestGNumerics::test_fft),
           "Test GFft");
    append(static_cast<pfunction>(&TestGNumerics::test_integral),
           "Test GIntegral");
    append(static_cast<pfunction>(&TestGNumerics::test_romberg_integration),
           "Test Romberg integration");
    append(static_cast<pfunction>(&TestGNumerics::test_adaptive_simpson_integration),
           "Test adaptive Simpson integration");
    append(static_cast<pfunction>(&TestGNumerics::test_gauss_kronrod_integration),
           "Test Gauss-Kronrod integration");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone test suite
 *
 * @return Pointer to deep copy of test suite.
 ***************************************************************************/
TestGNumerics* TestGNumerics::clone(void) const
{
    // Clone test suite
    return new TestGNumerics(*this);
}


/***********************************************************************//**
 * @brief Test mathematics functions
 ***************************************************************************/
void TestGNumerics::test_math(void)
{
    // Test acos function
    test_value(gammalib::acos(0.0), 1.57079633);
    test_value(gammalib::acos(0.37), 1.19178731);

    // Test atan2 function
    test_value(gammalib::atan2( 1.0, 2.0),  0.463647609001);
    test_value(gammalib::atan2( 1.0,-3.0),  2.819842099193);
    test_value(gammalib::atan2(-2.0,-3.0), -2.553590050042);
    test_value(gammalib::atan2(-2.0, 1.0), -1.107148717794);

    // Test cosd function
    test_value(gammalib::cosd(  0.0),  1.0);
    test_value(gammalib::cosd( 33.0),  0.83867057);
    test_value(gammalib::cosd(178.1), -0.99945022);

    // Test sind function
    test_value(gammalib::sind(  0.0), 0.0);
    test_value(gammalib::sind( 33.0), 0.54463904);
    test_value(gammalib::sind(178.1), 0.03315518);

    // Test tand function
    test_value(gammalib::tand(-10.0), -0.176326980708);
    test_value(gammalib::tand( 12.0),  0.21255656167);
    test_value(gammalib::tand(178.1), -0.0331734166041);

    // Test asind function
    test_value(gammalib::asind(0.0), 0.0);
    test_value(gammalib::asind(0.37), 21.7156172833);

    // Test acosd function
    test_value(gammalib::acosd(0.0), 90.0);
    test_value(gammalib::acosd(0.37), 68.28438272);

    // Test atan2d function
    test_value(gammalib::atan2d( 1.0, 2.0),   26.5650511771);
    test_value(gammalib::atan2d( 1.0,-3.0),  161.565051177);
    test_value(gammalib::atan2d(-2.0,-3.0), -146.309932474);
    test_value(gammalib::atan2d(-2.0, 1.0),  -63.4349488229);

    // Test sincosd function
    double s;
    double c;
    gammalib::sincosd(33.0, &s, &c);
    test_value(s, 0.54463904);
    test_value(c, 0.83867057);

    // Test logarithm of gamma function
    test_value(gammalib::gammln(0.5), 0.572364942925);
    test_value(gammalib::gammln(1.5), -0.120782237635);

    // Test error function
    test_value(gammalib::erf(0.0), 0.0);
    test_value(gammalib::erf(1.0), 0.8427008);
    test_value(gammalib::erf(2.0), 0.9953223);

    // Test complementary error function
    test_value(gammalib::erfc(0.0), 1.0);
    test_value(gammalib::erfc(0.1), 0.887537084);
    test_value(gammalib::erfc(0.8), 0.257899035);

    // Test inverse error function
    test_value(gammalib::erfinv(0.3), 0.2724627147);
    test_value(gammalib::erfinv(0.9), 1.1630871537);

    // Test modulo function
    test_value(gammalib::modulo(3.1, 2.8), 0.3);
    test_value(gammalib::modulo(-5.7, 1.37), 1.15);

    // Test power law integration
    test_value(gammalib::plaw_integral(2.0, 3.0, 9.0, 5.0), 29.112577);
    test_value(gammalib::plaw_integral(2.0, 9.0, 9.0, 2.0), 27.073393);
    test_value(gammalib::plaw_integral(2.0, 9.0, 5.0, 2.0), 12.471065);

    // Test Gaussian integration
    test_value(gammalib::gauss_integral(0.0, 1.0),  0.341344746069);
    test_value(gammalib::gauss_integral(-1.0, 1.0), 0.682689492137086);
    test_value(gammalib::gauss_integral(-2.0, 2.0), 0.954499736103642);
    test_value(gammalib::gauss_integral(-3.0, 3.0), 0.997300203936740);

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Test n-dimensional array
 ***************************************************************************/
void TestGNumerics::test_ndarray(void)
{
    // Test void constructor
    GNdarray nd0;
    test_value(nd0.size(), 0);
    test_value(nd0.dim(), 0);

    // Test 1-dimensional array
    GNdarray nd1(3);
    test_value(nd1.size(), 3);
    test_value(nd1.dim(), 1);
    test_value(nd1.shape()[0], 3);
    for (int i = 0; i < 3; ++i) {
        nd1(i) = double(i);
        test_value(nd1.at(i), double(i));
    }

    // Test 2-dimensional array
    GNdarray nd2(3,2);
    test_value(nd2.size(), 6);
    test_value(nd2.dim(), 2);
    test_value(nd2.shape()[0], 3);
    test_value(nd2.shape()[1], 2);
    for (int ix = 0; ix < 3; ++ix) {
        for (int iy = 0; iy < 2; ++iy) {
            double ref = double(ix)+double(iy);
            nd2(ix,iy) = ref;
            test_value(nd2.at(ix,iy), ref);
        }
    }
    // Test 3-dimensional array
    GNdarray nd3(2,2,2);
    test_value(nd3.size(), 8);
    test_value(nd3.dim(), 3);
    test_value(nd3.shape()[0], 2);
    test_value(nd3.shape()[1], 2);
    test_value(nd3.shape()[2], 2);
    for (int ix = 0; ix < 2; ++ix) {
        for (int iy = 0; iy < 2; ++iy) {
            for (int iz = 0; iz < 2; ++iz) {
                double ref    = double(ix)+double(iy)+double(iz);
                nd3(ix,iy,iz) = ref;
                test_value(nd3.at(ix,iy,iz), ref);
            }
        }
    }

    // Test n-dimensional array
    std::vector<int> n(5,2);
    GNdarray ndn(n);
    test_value(ndn.size(), 32);
    test_value(ndn.dim(), 5);
    for (int i = 0; i < 5; ++i) {
        test_value(ndn.shape()[i], 2);
    }
    std::vector<int> itest(5,1);
    ndn(itest) = 3.1415;
    test_value(ndn.at(itest), 3.1415);

    // Setup arrays for operator tests
    GNdarray nda(1,2);
    nda(0,0) = 3.0;
    nda(0,1) = 5.5;
    GNdarray ndb = nda;
    GNdarray ndc = nda + ndb;

    // Test operator=
    test_value(ndb(0,0), 3.0);
    test_value(ndb(0,1), 5.5);

    // Test operator== and operator!=
    test_assert(nda == ndb, "Test operator==");
    test_assert(nda != ndc, "Test operator!=");

    // Test array operator+=
    nda += ndb;
    test_value(nda(0,0),  6.0);
    test_value(nda(0,1), 11.0);

    // Test array operator-=
    nda -= ndb;
    test_value(nda(0,0), 3.0);
    test_value(nda(0,1), 5.5);

    // Test value operator+=
    nda += 1.5;
    test_value(nda(0,0), 4.5);
    test_value(nda(0,1), 7.0);

    // Test value operator-=
    nda -= 1.5;
    test_value(nda(0,0), 3.0);
    test_value(nda(0,1), 5.5);

    // Test value operator*=
    nda *= 3.0;
    test_value(nda(0,0),  9.0);
    test_value(nda(0,1), 16.5);

    // Test value operator/=
    nda /= 3.0;
    test_value(nda(0,0), 3.0);
    test_value(nda(0,1), 5.5);

    // Test operator-
    GNdarray ndd = -nda;
    test_value(ndd(0,0), -3.0);
    test_value(ndd(0,1), -5.5);

    // Test operator+(GNdarray,GNdarray)
    ndd = nda + ndb;
    test_value(ndd(0,0),  6.0);
    test_value(ndd(0,1), 11.0);

    // Test operator+(GNdarray,double) and operator+(double,GNdarray)
    ndd = nda + 1.5;
    test_value(ndd(0,0), 4.5);
    test_value(ndd(0,1), 7.0);
    ndd = 1.5 + nda;
    test_value(ndd(0,0), 4.5);
    test_value(ndd(0,1), 7.0);

    // Test operator-(GNdarray,GNdarray)
    ndd = nda - ndb;
    test_value(ndd(0,0), 0.0);
    test_value(ndd(0,1), 0.0);

    // Test operator-(GNdarray,double) and operator-(double,GNdarray)
    ndd = nda - 1.5;
    test_value(ndd(0,0), 1.5);
    test_value(ndd(0,1), 4.0);
    ndd = 1.5 - nda;
    test_value(ndd(0,0), -1.5);
    test_value(ndd(0,1), -4.0);

    // Test operator*(GNdarray,double) and operator*(double,GNdarray)
    ndd = nda * 3.0;
    test_value(ndd(0,0),  9.0);
    test_value(ndd(0,1), 16.5);
    ndd = 3.0 * nda;
    test_value(ndd(0,0),  9.0);
    test_value(ndd(0,1), 16.5);

    // Test operator/(GNdarray,double)
    ndd = nda / 3.0;
    test_value(ndd(0,0),     1.0);
    test_value(ndd(0,1), 5.5/3.0);

    // Test friend functions
    test_value(min(nda), 3.0);
    test_value(max(nda), 5.5);
    test_value(sum(nda), 8.5);
    test_value(acos(nda)(0,0), std::acos(3.0));
    test_value(acos(nda)(0,1), std::acos(5.5));
    test_value(acosh(nda)(0,0), acosh(3.0));
    test_value(acosh(nda)(0,1), acosh(5.5));
    test_value(asin(nda)(0,0), std::asin(3.0));
    test_value(asin(nda)(0,1), std::asin(5.5));
    test_value(asinh(nda)(0,0), asinh(3.0));
    test_value(asinh(nda)(0,1), asinh(5.5));
    test_value(atan(nda)(0,0), std::atan(3.0));
    test_value(atan(nda)(0,1), std::atan(5.5));
    test_value(atanh(nda)(0,0), atanh(3.0));
    test_value(atanh(nda)(0,1), atanh(5.5));
    test_value(cos(nda)(0,0), std::cos(3.0));
    test_value(cos(nda)(0,1), std::cos(5.5));
    test_value(cosh(nda)(0,0), std::cosh(3.0));
    test_value(cosh(nda)(0,1), std::cosh(5.5));
    test_value(exp(nda)(0,0), std::exp(3.0));
    test_value(exp(nda)(0,1), std::exp(5.5));
    test_value(abs(nda)(0,0), std::abs(3.0));
    test_value(abs(nda)(0,1), std::abs(5.5));
    test_value(log(nda)(0,0), std::log(3.0));
    test_value(log(nda)(0,1), std::log(5.5));
    test_value(log10(nda)(0,0), std::log10(3.0));
    test_value(log10(nda)(0,1), std::log10(5.5));
    test_value(sin(nda)(0,0), std::sin(3.0));
    test_value(sin(nda)(0,1), std::sin(5.5));
    test_value(sinh(nda)(0,0), std::sinh(3.0));
    test_value(sinh(nda)(0,1), std::sinh(5.5));
    test_value(sqrt(nda)(0,0), std::sqrt(3.0));
    test_value(sqrt(nda)(0,1), std::sqrt(5.5));
    test_value(tan(nda)(0,0), std::tan(3.0));
    test_value(tan(nda)(0,1), std::tan(5.5));
    test_value(tanh(nda)(0,0), std::tanh(3.0));
    test_value(tanh(nda)(0,1), std::tanh(5.5));
    test_value(pow(nda,3.0)(0,0), std::pow(3.0,3.0));
    test_value(pow(nda,3.0)(0,1), std::pow(5.5,3.0));

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Test Fast Fourier Transform
 ***************************************************************************/
void TestGNumerics::test_fft(void)
{
    // Test 1-dimensional FFT
    for (int n = 2; n < 12; ++n) {

        // Skip n=8,9,10 since they can be factorised
        if (n > 7 && n < 11) {
            continue;
        }

        // Allocate 1-dimensional array
        GNdarray array(n);

        // Set one elements
        array(n/2) = 1.0;

        // Allocate FFT
        GFft fft(array);

        // Backward transformation
        GNdarray back = fft.backward();

        // Check FFT results
        check_fft1(array, fft, back);

    } // endfor: looped over array lengths

    // Test 2-dimensional FFT
    GNdarray array2(8,9);
    array2(4,4) = 1.0;
    GFft     fft2(array2);
    GNdarray back2 = fft2.backward();
    check_fft2(array2, fft2, back2);

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Check 1-dimensional Fast Fourier Transform
 *
 * @param[in] array Input array.
 * @param[in] fft Fast Fourier Transform of input array.
 * @param[in] back Backward transform of FFT (should be equal to input array)
 ***************************************************************************/
void TestGNumerics::check_fft1(const GNdarray& array,
                               const GFft&     fft,
                               const GNdarray& back)
{
    // Get normalisation factor
    double norm = 1.0 / double(array.size());

    // Compute reference value and compare with FFT result
    for (int k = 0; k < array.size(); ++k) {
        double ref_real = 0.0;
        double ref_imag = 0.0;
        for (int n = 0; n < array.size(); ++n) {
            double x = -gammalib::twopi * k * n * norm;
            ref_real += array(n) * std::cos(x);
            ref_imag += array(n) * std::sin(x);
        }
        std::complex<double> ref(ref_real, ref_imag);
        test_value(fft(k), ref);
    }

    // Test backward transformation
    for (int k = 0; k < array.size(); ++k) {
        test_value(back(k), array(k));
    }

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Check 2-dimensional Fast Fourier Transform
 *
 * @param[in] array Input array.
 * @param[in] fft Fast Fourier Transform of input array.
 * @param[in] back Backward transform of FFT (should be equal to input array)
 ***************************************************************************/
void TestGNumerics::check_fft2(const GNdarray& array,
                               const GFft&     fft,
                               const GNdarray& back)
{
    // Get normalisation factor
    double normM = 1.0 / double(array.shape()[0]);
    double normN = 1.0 / double(array.shape()[1]);

    // Compute reference value and compare with FFT result
    for (int k = 0; k < array.shape()[0]; ++k) {
        for (int l = 0; l < array.shape()[1]; ++l) {
            double ref_real = 0.0;
            double ref_imag = 0.0;
            for (int n = 0; n < array.shape()[1]; ++n) {
                double sum_real = 0.0;
                double sum_imag = 0.0;
                for (int m = 0; m < array.shape()[0]; ++m) {
                    double x = -gammalib::twopi * m * k * normM;
                    sum_real += array(m,n) * std::cos(x);
                    sum_imag += array(m,n) * std::sin(x);
                }
                double y = -gammalib::twopi * n * l * normN;
                ref_real += (sum_real * std::cos(y) - sum_imag * std::sin(y));
                ref_imag += (sum_real * std::sin(y) + sum_imag * std::cos(y));
            }
            std::complex<double> ref(ref_real, ref_imag);
            test_value(fft(k,l), ref);
        }
    }

    // Test backward transformation
    for (int k = 0; k < array.size(); ++k) {
        test_value(back(k), array(k));
    }

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Test model parameter handling
 ***************************************************************************/
void TestGNumerics::test_integral(void)
{
    // Test integral and integrand allocation
    Gauss     integrand(m_sigma);
    GIntegral integral(&integrand);

    // Exit test
    return;
}

/***********************************************************************//**
 * @brief Test Romberg integration
 ***************************************************************************/
void TestGNumerics::test_romberg_integration(void)
{
    // Set-up integral
    Gauss     integrand(m_sigma);
    GIntegral integral(&integrand);

    // Integrate over the entire Gaussian
    double result = integral.romberg(-10.0*m_sigma, 10.0*m_sigma);
    test_value(result,1.0,1.0e-6,"","Gaussian integral is not 1.0 (integral="+gammalib::str(result)+")");

    // Test [-1sigma, 1sigma]
    result = integral.romberg(-m_sigma, m_sigma);
    test_value(result,0.68268948130801355,1.0e-6,"","Gaussian integral is not 0.682689 (difference="+gammalib::str((result-0.68268948130801355))+")");

    // Test [0.0, 1sigma]
    result = integral.romberg(0.0, m_sigma);
    test_value(result,0.3413447460687748,1.0e-6,"","Gaussian integral is not 0.341345 (difference="+gammalib::str((result-0.3413447460687748))+")");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test adaptive Simpson integration
 ***************************************************************************/
void TestGNumerics::test_adaptive_simpson_integration(void)
{
    // Set-up integral
    Gauss     integrand(m_sigma);
    GIntegral integral(&integrand);

    // Integrate over the entire Gaussian
    double result = integral.adaptive_simpson(-10.0*m_sigma, 10.0*m_sigma);
/*
std::cout << result << std::endl;
std::cout << result-1.0 << std::endl;
std::cout << integral << std::endl;
*/
    test_value(result,1.0,1.0e-6,"","Gaussian integral is not 1.0 (integral="+gammalib::str(result)+")");

    // Test [-1sigma, 1sigma]
    result = integral.adaptive_simpson(-m_sigma, m_sigma);
/*
std::cout << result << std::endl;
std::cout << result-0.68268948130801355 << std::endl;
std::cout << integral << std::endl;
*/
    test_value(result,0.68268948130801355,1.0e-6,"","Gaussian integral is not 0.682689 (difference="+gammalib::str((result-0.68268948130801355))+")");

    // Test [0.0, 1sigma]
    result = integral.adaptive_simpson(0.0, m_sigma);
/*
std::cout << result << std::endl;
std::cout << result-0.3413447460687748 << std::endl;
std::cout << integral << std::endl;
*/
    test_value(result,0.3413447460687748,1.0e-6,"","Gaussian integral is not 0.341345 (difference="+gammalib::str((result-0.3413447460687748))+")");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test Gauss-Kronrod integration
 ***************************************************************************/
void TestGNumerics::test_gauss_kronrod_integration(void)
{
    // Set-up integral
    Gauss     integrand(m_sigma);
    GIntegral integral(&integrand);

    // Integrate over the entire Gaussian
    double result = integral.gauss_kronrod(-10.0*m_sigma, 10.0*m_sigma);
/*
std::cout << result << std::endl;
std::cout << result-1.0 << std::endl;
std::cout << integral << std::endl;
*/
    test_value(result,1.0,1.0e-6,"","Gaussian integral is not 1.0 (integral="+gammalib::str(result)+")");

    // Test [-1sigma, 1sigma]
    result = integral.gauss_kronrod(-m_sigma, m_sigma);
/*
std::cout << result << std::endl;
std::cout << result-0.68268948130801355 << std::endl;
std::cout << integral << std::endl;
*/
    test_value(result,0.68268948130801355,1.0e-6,"","Gaussian integral is not 0.682689 (difference="+gammalib::str((result-0.68268948130801355))+")");

    // Test [0.0, 1sigma]
    result = integral.gauss_kronrod(0.0, m_sigma);
/*
std::cout << result << std::endl;
std::cout << result-0.3413447460687748 << std::endl;
std::cout << integral << std::endl;
*/
    test_value(result,0.3413447460687748,1.0e-6,"","Gaussian integral is not 0.341345 (difference="+gammalib::str((result-0.3413447460687748))+")");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Main test function
 ***************************************************************************/
int main(void)
{
    GTestSuites testsuite("GNumerics");

    bool was_successful=true;

    //Create a test suite
    TestGNumerics test;

    //Append to the container
    testsuite.append(test);

    //Run
    was_successful=testsuite.run();

    //save xml report
    testsuite.save("reports/GNumerics.xml");

    // Return
    return was_successful ? 0:1;
}
