/***************************************************************************
 *                test_GModel.hpp  -  test model class                     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2013 by Jean-Baptiste Cayrou                        *
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
 * @file test_GModel.cpp
 * @brief Test model classes
 * @author Jean-Baptiste Cayrou
 */

#ifndef TEST_GMODEL_HPP
#define TEST_GMODEL_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"


/***********************************************************************//**
 * @class TestGModel
 *
 * @brief Test suite for model module testing
 ***************************************************************************/
class TestGModel : public GTestSuite
{
public:
    // Constructors and destructors
    TestGModel(void) : GTestSuite() {}
    virtual ~TestGModel(void) {}

    // Methods
    virtual void set(void);
    void    test_model_par(void);
    void    test_sky_model(void);
    void    test_point_source(void);
    void    test_diffuse_const(void);
    void    test_diffuse_cube(void);
    void    test_diffuse_map(void);
    void    test_radial_disk(void);
    void    test_radial_gauss(void);
    void    test_radial_shell(void);
    void    test_elliptical_disk(void);
    void    test_spatial_model(void);
    void    test_const(void);
    void    test_plaw(void);
    void    test_plaw2(void);
    void    test_eplaw(void);
    void    test_logparabola(void);
    void    test_nodes(void);
    void    test_filefct(void);
    void    test_spectral_model(void);
    void    test_temp_const(void);
    void    test_models(void);
    void    test_model_registry(void);

private:        
    // Private methods
    void test_xml_model(const std::string& name, const std::string& filename);
    
    // Private attributes
    std::string m_map_file;
    std::string m_cube_file;
    std::string m_xml_file;
    std::string m_xml_model_point_const;
    std::string m_xml_model_point_plaw;
    std::string m_xml_model_point_plaw2;
    std::string m_xml_model_point_eplaw;
    std::string m_xml_model_point_logparabola;
    std::string m_xml_model_point_nodes;
    std::string m_xml_model_point_filefct;
    std::string m_xml_model_diffuse_const;
    std::string m_xml_model_diffuse_cube;
    std::string m_xml_model_diffuse_map;
    std::string m_xml_model_radial_disk;
    std::string m_xml_model_radial_gauss;
    std::string m_xml_model_radial_shell;
    std::string m_xml_model_elliptical_disk;
};

#endif /* TEST_GMODEL_HPP */
