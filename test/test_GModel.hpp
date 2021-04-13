/***************************************************************************
 *                test_GModel.hpp - Test model module                      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2021 by Jean-Baptiste Cayrou                        *
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
    virtual void        set(void);
    virtual TestGModel* clone(void) const;
    virtual std::string classname(void) const { return "TestGModel"; }
    void                test_model_par(void);
    void                test_model_association(void);
    void                test_model_associations(void);
    void                test_sky_model(void);
    void                test_point_source(void);
    void                test_diffuse_const(void);
    void                test_diffuse_cube(void);
    void                test_diffuse_map(void);
    void                test_radial_disk(void);
    void                test_radial_ring(void);
    void                test_radial_gauss(void);
    void                test_radial_shell(void);
    void                test_elliptical_disk(void);
    void                test_elliptical_gauss(void);
    void                test_spatial_composite(void);
    void                test_spatial_model(void);
    void                test_const(void);
    void                test_gauss(void);
    void                test_plaw(void);
    void                test_plaw_phflux(void);
    void                test_plaw_eflux(void);
    void                test_eplaw(void);
    void                test_einvplaw(void);
    void                test_exponential(void);
    void                test_bplaw(void);
    void                test_smoothbplaw(void);
    void                test_supeplaw(void);
    void                test_logparabola(void);
    void                test_multiplicative(void);
    void                test_bins(void);
    void                test_nodes(void);
    void                test_filefct(void);
    void                test_spectral_composite(void);
    void                test_spectral_model(void);
    void                test_table(void);
    void                test_temp_const(void);
    void                test_temp_lightcurve(void);
    void                test_temp_phasecurve(void);
    void                test_model(void);
    void                test_models(void);
    void                test_model_registry(void);
    void                test_legacy_model_radial_disk(void);
    void                test_legacy_model_radial_gauss(void);
    void                test_legacy_model_radial_shell(void);
    void                test_legacy_model_elliptical_gauss(void);
    void                test_legacy_model_diffuse_const(void);
    void                test_legacy_model_diffuse_map(void);
    void                test_legacy_model_diffuse_cube(void);
    void                test_legacy_model_point_const(void);
    void                test_legacy_model_point_plaw(void);
    void                test_legacy_model_point_plaw2(void);
    void                test_legacy_model_point_smoothbplaw(void);
    void                test_legacy_model_point_eplaw(void);
    void                test_legacy_model_point_supeplaw(void);
    void                test_legacy_model_point_logparabola(void);
    void                test_flux(void);

private:
    // Private methods
    void test_xml_model(const std::string& name, const std::string& filename);
    void test_sky_model_content(const GModelSky& model);

    // Test files
    std::string m_map_file;
    std::string m_cube_file;
    std::string m_filefct;
    std::string m_table;
    std::string m_temp_lightcurve;
    std::string m_temp_phasecurve;
    std::string m_xml_file;
    std::string m_assoc_file;

    // Model definiton XML files
    std::string m_xml_model_point_const;
    std::string m_xml_model_point_gauss;
    std::string m_xml_model_point_plaw;
    std::string m_xml_model_point_plaw_phflux;
    std::string m_xml_model_point_plaw_eflux;
    std::string m_xml_model_point_eplaw;
    std::string m_xml_model_point_einvplaw;
    std::string m_xml_model_point_bplaw;
    std::string m_xml_model_point_smoothbplaw;
    std::string m_xml_model_point_supeplaw;
    std::string m_xml_model_point_logparabola;
    std::string m_xml_model_point_bins;
    std::string m_xml_model_point_nodes;
    std::string m_xml_model_point_filefct;
    std::string m_xml_model_point_table;
    std::string m_xml_point_multiplicative;
    std::string m_xml_point_exponential;
    std::string m_xml_model_spectral_composite;
    std::string m_xml_model_diffuse_const;
    std::string m_xml_model_diffuse_cube;
    std::string m_xml_model_diffuse_map;
    std::string m_xml_model_radial_disk;
    std::string m_xml_model_radial_ring;
    std::string m_xml_model_radial_gauss;
    std::string m_xml_model_radial_shell;
    std::string m_xml_model_elliptical_disk;
    std::string m_xml_model_elliptical_gauss;
    std::string m_xml_model_spatial_composite;
    std::string m_xml_model_point_temp_lightcurve;
    std::string m_xml_model_point_temp_phasecurve;

    // Legacy model definition XML files
    std::string m_xml_legacy_radial_disk;
    std::string m_xml_legacy_radial_gauss;
    std::string m_xml_legacy_radial_shell;
    std::string m_xml_legacy_elliptical_gauss;
    std::string m_xml_legacy_diffuse_const;
    std::string m_xml_legacy_diffuse_map;
    std::string m_xml_legacy_diffuse_cube;
    std::string m_xml_legacy_point_const;
    std::string m_xml_legacy_point_plaw;
    std::string m_xml_legacy_point_plaw2;
    std::string m_xml_legacy_point_smoothbplaw;
    std::string m_xml_legacy_point_eplaw;
    std::string m_xml_legacy_point_supeplaw;
    std::string m_xml_legacy_point_logparabola;
};

#endif /* TEST_GMODEL_HPP */
