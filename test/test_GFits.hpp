/***************************************************************************
 *              test_GFits.hpp  -   test fits class                        *
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

#ifndef TEST_GFITS_HPP
#define TEST_GFITS_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"
#include <iostream>                           // cout, cerr

class TestGFits : public GTestSuite
{
    public:
        // Constructors and destructors
        TestGFits(void) : GTestSuite(){ return; }
        virtual ~TestGFits(void){ return; }

        // Methods
        virtual void set(void);
        void         test_create(void);
        void         test_image_byte(void);
        void         test_image_ushort(void);
        void         test_image_short(void);
        void         test_image_ulong(void);
        void         test_image_long(void);
        void         test_image_longlong(void);
        void         test_image_float(void);
        void         test_image_double(void);
        void         test_bintable_bit(void);
        void         test_bintable_logical(void);
        void         test_bintable_string(void);
        void         test_bintable_double(void);
        void         test_bintable_float(void);
        void         test_bintable_ushort(void);
        void         test_bintable_short(void);
        void         test_bintable_ulong(void);
        void         test_bintable_long(void);
        void         test_bintable_longlong(void);

    // Private methods
    private:
        int equal(double val, double ref, double eps);
        int fequal(double val, double ref);
        int dequal(double val, double ref);


};

#endif /* TEST_GFITS_HPP */
