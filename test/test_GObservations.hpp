/***************************************************************************
 *              test_GObservation.hpp  -  test GObservation class          *
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

#ifndef TEST_GOBSERVATION_HPP
#define TEST_GOBSERVATION_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"
#include <iostream>                           // cout, cerr
#include "testinst/GTestLib.hpp"              // Test instrument

#ifdef _OPENMP
class TestGObservation : public GTestSuite
{
    public:
        // Constructors and destructors
        TestGObservation(const std::string& name);
        virtual ~TestGObservation(void);

        // Methods
        GModelPar& test_observations_optimizer(int mode=0);
        void test_observations_optimizer_unbinned_1();
        void test_observations_optimizer_unbinned_10();
        void test_observations_optimizer_binned_1();
        void test_observations_optimizer_binned_10();
        virtual void set(void);

    // Private members
    private:
};
#endif

#endif /* TEST_GOPTIMIZER_HPP */
