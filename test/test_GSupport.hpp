/***************************************************************************
 *              test_GSupport.hpp  -   test support class                  *
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

#ifndef TEST_GSUPPORT_HPP
#define TEST_GSUPPORT_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"
#include <iostream>                           // cout, cerr

class TestGSupport : public GTestSuite
{
    public:
        // Constructors and destructors
        TestGSupport(void) : GTestSuite(){ return; }
        virtual ~TestGSupport(void){ return; }

        // Methods
        virtual void set(void);
        void test_expand_env(void);

    // Private members
    private:
};

#endif /* TEST_GSUPPORT_HPP */
