/***************************************************************************
 *                test_GModel.hpp  -  test model class                     *
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

#ifndef TEST_GMODEL_HPP
#define TEST_GMODEL_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"
#include "GTools.hpp"

class TestGModel : public GTestSuite
{
    public:
        // Constructors and destructors
        TestGModel(void) : GTestSuite(){ return; }
        virtual ~TestGModel(void){ return; }

        // Methods
        virtual void set(void);

        void test_model_par(void);
        void test_model(void);
        void test_models(void);
        void test_spectral_model(void);
        void test_spacial_model(void);
    // Private attributes
    private:
        std::string m_xml_file;
        std::string m_xml_model_point_nodes;
        std::string m_xml_model_spatial_map;
        
    //private methods
        void test_xml_model(const std::string &name, const std::string &filename);

};

#endif /* TEST_GMODEL_HPP */
