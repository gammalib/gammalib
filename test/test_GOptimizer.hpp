/***************************************************************************
 *              test_GOptimizer.hpp  -  test GOptimizer class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2012 by Jurgen Knodlseder                           *
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

#ifndef TEST_GOPTIMIZER_HPP
#define TEST_GOPTIMIZER_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"
#include <iostream>                           // cout, cerr
#include <stdexcept>                          // std::exception

class TestGOptimizer : public GTestSuite
{
    public:
        // Constructors and destructors
        TestGOptimizer(void) : GTestSuite(){ return; }
        virtual ~TestGOptimizer(void){ return; }

        // Methods
        virtual void set(void);
        void test_optimizer(void);

    // Private members
    private:
};


#endif /* TEST_GOPTIMIZER_HPP */
