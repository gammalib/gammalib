/***************************************************************************
 *                   test_GModel.cpp - test GModel class                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2016 by Juergen Knoedlseder                         *
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
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iostream>
#include <ostream>
#include <stdexcept>
#include <cstdlib>     // getenv
#include "test_GModel.hpp"
#include "GTools.hpp"


/***********************************************************************//**
 * @brief Set tests
 ***************************************************************************/
void TestGModel::set(void)
{
    // Set name
    name("GModel");

    // Set test data directory
    std::string datadir = std::getenv("TEST_DATA");

    // Set test files
    m_map_file  = datadir + "/cena_lobes_parkes.fits";
    m_cube_file = datadir + "/test_cube.fits";
    m_filefct   = datadir + "/filefunction.txt";
    m_xml_file  = datadir + "/crab.xml";

    // Set model definiton XML files
    m_xml_model_point_const        = datadir + "/model_point_const.xml";
    m_xml_model_point_gauss        = datadir + "/model_point_gauss.xml";
    m_xml_model_point_plaw         = datadir + "/model_point_plaw.xml";
    m_xml_model_point_plaw_phflux  = datadir + "/model_point_plaw_phflux.xml";
    m_xml_model_point_plaw_eflux   = datadir + "/model_point_plaw_eflux.xml";
    m_xml_model_point_eplaw        = datadir + "/model_point_eplaw.xml";
    m_xml_model_point_einvplaw     = datadir + "/model_point_einvplaw.xml";
    m_xml_model_point_bplaw        = datadir + "/model_point_bplaw.xml";
    m_xml_model_point_supeplaw     = datadir + "/model_point_supeplaw.xml";
    m_xml_model_point_logparabola  = datadir + "/model_point_logparabola.xml";
    m_xml_model_point_nodes        = datadir + "/model_point_nodes.xml";
    m_xml_model_point_filefct      = datadir + "/model_point_filefct.xml";
    m_xml_model_diffuse_const      = datadir + "/model_diffuse_const.xml";
    m_xml_model_diffuse_cube       = datadir + "/model_diffuse_cube.xml";
    m_xml_model_diffuse_map        = datadir + "/model_diffuse_map.xml";
    m_xml_model_radial_disk        = datadir + "/model_radial_disk.xml";
    m_xml_model_radial_gauss       = datadir + "/model_radial_gauss.xml";
    m_xml_model_radial_shell       = datadir + "/model_radial_shell.xml";
    m_xml_model_elliptical_disk    = datadir + "/model_elliptical_disk.xml";
    m_xml_model_elliptical_gauss   = datadir + "/model_elliptical_gauss.xml";

    // Set legacy model definition XML files
    m_xml_legacy_radial_disk       = datadir + "/legacy_radial_disk.xml";
    m_xml_legacy_radial_gauss      = datadir + "/legacy_radial_gauss.xml";
    m_xml_legacy_radial_shell      = datadir + "/legacy_radial_shell.xml";
    m_xml_legacy_elliptical_gauss  = datadir + "/legacy_elliptical_gauss.xml";
    m_xml_legacy_diffuse_const     = datadir + "/legacy_diffuse_const.xml";
    m_xml_legacy_diffuse_map       = datadir + "/legacy_diffuse_map.xml";
    m_xml_legacy_diffuse_cube      = datadir + "/legacy_diffuse_cube.xml";
    m_xml_legacy_point_const       = datadir + "/legacy_point_const.xml";
    m_xml_legacy_point_plaw        = datadir + "/legacy_point_plaw.xml";
    m_xml_legacy_point_plaw2       = datadir + "/legacy_point_plaw2.xml";
    m_xml_legacy_point_eplaw       = datadir + "/legacy_point_eplaw.xml";
    m_xml_legacy_point_supeplaw    = datadir + "/legacy_point_supeplaw.xml";
    m_xml_legacy_point_logparabola = datadir + "/legacy_point_logparabola.xml";

    // Append tests
    append(static_cast<pfunction>(&TestGModel::test_model_par),
           "Test GModelPar");
    append(static_cast<pfunction>(&TestGModel::test_sky_model),
           "Test GModelSky");

    // Append spatial model tests
    append(static_cast<pfunction>(&TestGModel::test_point_source),
           "Test GModelSpatialPointSource");
    append(static_cast<pfunction>(&TestGModel::test_radial_disk),
           "Test GModelSpatialRadialDisk");
    append(static_cast<pfunction>(&TestGModel::test_radial_gauss),
           "Test GModelSpatialRadialGauss");
    append(static_cast<pfunction>(&TestGModel::test_radial_shell),
           "Test GModelSpatialRadialShell");
    append(static_cast<pfunction>(&TestGModel::test_elliptical_disk),
           "Test GModelSpatialEllipticalDisk");
    append(static_cast<pfunction>(&TestGModel::test_elliptical_gauss),
           "Test GModelSpatialEllipticalGauss");
    append(static_cast<pfunction>(&TestGModel::test_diffuse_const),
           "Test GModelSpatialDiffuseConst");
    append(static_cast<pfunction>(&TestGModel::test_diffuse_cube),
           "Test GModelSpatialDiffuseCube");
    append(static_cast<pfunction>(&TestGModel::test_diffuse_map),
           "Test GModelSpatialDiffuseMap");
    append(static_cast<pfunction>(&TestGModel::test_spatial_model),
           "Test spatial model XML I/O");

    // Append spectral model tests
    append(static_cast<pfunction>(&TestGModel::test_const),
           "Test GModelSpectralConst");
    append(static_cast<pfunction>(&TestGModel::test_gauss),
           "Test GModelSpectralGauss");
    append(static_cast<pfunction>(&TestGModel::test_plaw),
           "Test GModelSpectralPlaw");
    append(static_cast<pfunction>(&TestGModel::test_plaw_phflux),
           "Test GModelSpectralPlawPhotonFkux");
    append(static_cast<pfunction>(&TestGModel::test_plaw_eflux),
           "Test GModelSpectralPlawEnergyFlux");
    append(static_cast<pfunction>(&TestGModel::test_eplaw),
           "Test GModelSpectralExpPlaw");
    append(static_cast<pfunction>(&TestGModel::test_einvplaw),
           "Test GModelSpectralExpInvPlaw");
    append(static_cast<pfunction>(&TestGModel::test_supeplaw),
           "Test GModelSpectralSuperExpPlaw");
    append(static_cast<pfunction>(&TestGModel::test_bplaw),
           "Test GModelSpectralBrokenPlaw");
    append(static_cast<pfunction>(&TestGModel::test_logparabola),
           "Test GModelSpectralLogParabola");
    append(static_cast<pfunction>(&TestGModel::test_nodes),
           "Test GModelSpectralNodes");
    append(static_cast<pfunction>(&TestGModel::test_filefct),
           "Test GModelSpectralFunc");
    append(static_cast<pfunction>(&TestGModel::test_spectral_model),
           "Test spectral model XML I/O");

    // Append temporal model tests
    append(static_cast<pfunction>(&TestGModel::test_temp_const),
           "Test GModelTemporalConst");

    // Append genereric model tests
    append(static_cast<pfunction>(&TestGModel::test_model), "Test GModel");

    // Append model container tests
    append(static_cast<pfunction>(&TestGModel::test_models), "Test GModels");

    // Append model registry tests
    append(static_cast<pfunction>(&TestGModel::test_model_registry),
           "Test model registries");

    // Append legacy model tests
    #if defined(G_LEGACY_XML_FORMAT)
    append(static_cast<pfunction>(&TestGModel::test_legacy_model_radial_disk),
           "Test GModelSpatialRadialDisk legacy model");
    append(static_cast<pfunction>(&TestGModel::test_legacy_model_radial_gauss),
           "Test GModelSpatialRadialGauss legacy model");
    append(static_cast<pfunction>(&TestGModel::test_legacy_model_radial_shell),
           "Test GModelSpatialRadialShell legacy model");
    append(static_cast<pfunction>(&TestGModel::test_legacy_model_elliptical_gauss),
           "Test GModelSpatialEllipticalGauss legacy model");
    append(static_cast<pfunction>(&TestGModel::test_legacy_model_diffuse_const),
           "Test GModelSpatialDiffuseConst legacy model");
    append(static_cast<pfunction>(&TestGModel::test_legacy_model_diffuse_map),
           "Test GModelSpatialDiffuseMap legacy model");
    append(static_cast<pfunction>(&TestGModel::test_legacy_model_diffuse_cube),
           "Test GModelSpatialDiffuseCube legacy model");
    append(static_cast<pfunction>(&TestGModel::test_legacy_model_point_const),
           "Test GModelSpectralConst legacy model");
    append(static_cast<pfunction>(&TestGModel::test_legacy_model_point_plaw),
           "Test GModelSpectralPlaw legacy model");
    append(static_cast<pfunction>(&TestGModel::test_legacy_model_point_plaw2),
           "Test GModelSpectralPlaw2 legacy model");
    append(static_cast<pfunction>(&TestGModel::test_legacy_model_point_eplaw),
           "Test GModelSpectralExpPlaw legacy model");
    append(static_cast<pfunction>(&TestGModel::test_legacy_model_point_supeplaw),
           "Test GModelSpectralSuperExpPlaw legacy model");
    append(static_cast<pfunction>(&TestGModel::test_legacy_model_point_logparabola),
           "Test GModelSpectralLogParabola legacy model");
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone test suite
 *
 * @return Pointer to deep copy of test suite.
 ***************************************************************************/
TestGModel* TestGModel::clone(void) const
{
    // Clone test suite
    return new TestGModel(*this);
}


/***********************************************************************//**
 * @brief Test GModelPar class
 ***************************************************************************/
void TestGModel::test_model_par(void)
{
    // Test void constructor
    test_try("Test void constructor");
    try {
        GModelPar par;
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }
    
    // Test parameter constructor (value version)
    test_try("Test parameter constructor");
    try {
        GModelPar par("Test parameter", 47.0);
        test_value(par.value(), 47.0);
        test_value(par.error(), 0.0);
        test_value(par.gradient(), 0.0);
        test_value(par.min(), 0.0);
        test_value(par.max(), 0.0);
        test_assert(!par.has_min(), "Parameter shall have no minimum.");
        test_assert(!par.has_max(), "Parameter shall have no maximum.");
        test_assert(!par.has_range(), "Parameter shall have no range.");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }
    
    // Test invalid parameter constructor (factor & scale version)
    test_try("Test invalid parameter constructor");
    try {
        GModelPar par("Test parameter", 1.0, 0.0);
        test_try_failure("Parameter constructor with zero scale factor"
                         " shall throw an exception.");
    }
    catch (GException::invalid_argument &e) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test valid parameter constructor (factor & scale version)
    test_try("Test invalid parameter constructor");
    try {
        GModelPar par("Test parameter", 47.01, 2.0);
        par.factor_error(2.003);
        par.factor_gradient(51.0);
        par.unit("MeV");
        par.free();
        test_value(par.value(), 94.02);
        test_value(par.error(), 4.006);
        test_value(par.gradient(), 25.5);
        test_value(par.factor_value(), 47.01);
        test_value(par.factor_error(), 2.003);
        test_value(par.factor_gradient(), 51.0);
        test_value(par.scale(), 2.0);
        test_assert(par.name() == "Test parameter", "Parameter name");
        test_assert(par.unit() == "MeV", "Parameter unit");
        test_assert(par.is_free(), "Parameter freezing");
        test_assert(par.print() == 
                    "  Test parameter ...........: 94.02 +/- 4.006 MeV"
                    " (free,scale=2)", "Parameter printing");
        GModelPar par2("Another test parameter", 3.14, 3.0);
        test_value(par2.value(), 9.42);
        test_value(par2.factor_value(), 3.14);
        test_value(par2.scale(), 3.0);
        test_assert(par2.name() == "Another test parameter", "Parameter name");
        test_assert(par2.is_free(), "Parameter freezing");
        test_assert(par2.print() == 
                    "  Another test parameter ...: 9.42 +/- 0  (free,scale=3)",
                    "Parameter printing");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test boundary handling 1
    test_try("Test boundary handling (1/4)");
    try {
        GModelPar par("Test boundary", 1.0);
        test_assert(!par.has_min(), "Parameter shall have no minimum.");
        test_assert(!par.has_max(), "Parameter shall have no maximum.");
        test_assert(!par.has_range(), "Parameter shall have no range.");
        par.min(0.5);
        test_assert(par.has_min(), "Parameter shall have minimum.");
        test_assert(!par.has_max(), "Parameter shall have no maximum.");
        test_assert(!par.has_range(), "Parameter shall have no range.");
        par.max(2.0);
        test_assert(par.has_min(), "Parameter shall have minimum.");
        test_assert(par.has_max(), "Parameter shall have maximum.");
        test_assert(par.has_range(), "Parameter shall have range.");
        par.value(5.0);
        test_try_failure("Setting a parameter outside boundaries shall"
                         " generate an exception.");
    }
    catch (GException::invalid_argument &e) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test boundary handling 2/4
    test_try("Test boundary handling (2/4)");
    try {
        GModelPar par("Test boundary", 1.0);
        par.min(2.0);
        test_try_failure("Setting the minimum boundary that is larger"
                         " than the parameter value shall generate an"
                         " exception.");
    }
    catch (GException::invalid_argument &e) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test boundary handling 3/4
    test_try("Test boundary handling (3/4)");
    try {
        GModelPar par("Test boundary", 1.0);
        par.max(0.5);
        test_try_failure("Setting the maximum boundary that is smaller"
                         " than the parameter value shall generate an"
                         " exception.");
    }
    catch (GException::invalid_argument &e) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test boundary handling 4/4
    test_try("Test boundary handling (4/4)");
    try {
        GModelPar par("Test boundary", 1.0);
        par.range(2.0, 0.5);
        test_try_failure("Setting the minimum boundary that is larger"
                         " than the maximum boundary shall generate an"
                         " exception.");
    }
    catch (GException::invalid_argument &e) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test autoscaling
    GModelPar par("Test parameter", 3.0);
    par.scale(1.0);
    par.error(3.0);
    par.gradient(3.0);
    par.min(3.0);
    par.max(3.0);
    test_value(par.scale(), 1.0);
    test_value(par.factor_value(), 3.0);
    test_value(par.factor_error(), 3.0);
    test_value(par.factor_gradient(), 3.0);
    test_value(par.factor_min(), 3.0);
    test_value(par.factor_max(), 3.0);
    par.autoscale();
    test_value(par.scale(), 3.0);
    test_value(par.factor_value(), 1.0);
    test_value(par.factor_error(), 1.0);
    test_value(par.factor_gradient(), 9.0);
    test_value(par.factor_min(), 1.0);
    test_value(par.factor_max(), 1.0);
    //
    par.clear();
    par.value(5.0e-16);
    par.error(1.0e-16);
    par.gradient(1.0e-16);
    par.min(3.0e-16);
    par.max(8.0e-16);
    par.autoscale();
    test_value(par.scale(), 5.0e-16);
    test_value(par.value(), 5.0e-16);
    test_value(par.error(), 1.0e-16);
    test_value(par.gradient(), 1.0e-16);
    test_value(par.min(), 3.0e-16);
    test_value(par.max(), 8.0e-16);
    //
    par.clear();
    par.value(-5.0e-16);
    par.error(1.0e-16);
    par.gradient(2.0e-16);
    par.min(-8.0e-16);
    par.max(-3.0e-16);
    par.autoscale();
    test_value(par.scale(), -5.0e-16);
    test_value(par.value(), -5.0e-16);
    test_value(par.error(), 1.0e-16);
    test_value(par.gradient(), 2.0e-16);
    test_value(par.min(), -3.0e-16);     // Inverted due to autoscale ???
    test_value(par.max(), -8.0e-16);     // Inverted due to autoscale ???

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Test GModelSky class
 ***************************************************************************/
void TestGModel::test_sky_model(void)
{
    // Test void constructor
    GModelSky sky1;
    test_value(sky1.size(), 0, "Void constructor has no parameters");
    test_assert(sky1.spatial() == NULL,
                "Void constructor has no spatial component");
    test_assert(sky1.spectral() == NULL,
                "Void constructor has no spectral component");
    test_assert(sky1.temporal() == NULL,
                "Void constructor has no temporal component");

    // Test type constructor
    GModelSky sky2("My type");
    test_value(sky2.size(), 0, "Type constructor has no parameters");
    test_assert(sky2.spatial() == NULL,
                "Type constructor has no spatial component");
    test_assert(sky2.spectral() == NULL,
                "Type constructor has no spectral component");
    test_assert(sky2.temporal() == NULL,
                "Type constructor has no temporal component");

    // Setup components for constructor tests
    GModelSpatialPointSource  spat_ptsrc(10.0, -9.0);
    GModelSpatialDiffuseConst spat_const(3.0);
    GModelSpectralPlaw        spec_plaw(1.0, -1.0, GEnergy(1.0, "TeV"));
    GModelSpectralConst       spec_const(2.0);
    GModelTemporalConst       temp_const(4.0);

    // Test spatial and spectral component constructor
    GModelSky sky3(spat_ptsrc, spec_const);
    test_value(sky3.size(), 4); // 2 spatial, 1 spectral, 1 temporal
    test_value(sky3.spatial()->classname(),  "GModelSpatialPointSource");
    test_value(sky3.spectral()->classname(), "GModelSpectralConst");
    test_value(sky3.temporal()->classname(), "GModelTemporalConst");

    // Test component constructor
    GModelSky sky4(spat_ptsrc, spec_const, temp_const);
    test_value(sky4.size(), 4); // 2 spatial, 1 spectral, 1 temporal
    test_value(sky4.spatial()->classname(),  "GModelSpatialPointSource");
    test_value(sky4.spectral()->classname(), "GModelSpectralConst");
    test_value(sky4.temporal()->classname(), "GModelTemporalConst");

    // Test spatial method
    sky4.spatial(&spat_const);
    test_value(sky4.size(), 3); // 1 spatial, 1 spectral, 1 temporal
    test_value(sky4.spatial()->classname(), "GModelSpatialDiffuseConst");

    // Test spectral method
    sky4.spectral(&spec_plaw);
    test_value(sky4.size(), 5); // 1 spatial, 3 spectral, 1 temporal
    test_value(sky4.spectral()->classname(), "GModelSpectralPlaw");
    
    // Test XML constructor
    GXml         xml1(m_xml_file);
    GXmlElement* element1 = xml1.element(0)->element(0);
    GModelSky    sky5(*element1);
    test_value(sky5.size(), 6, "Check number of model parameters");
    test_value(sky5.name(), "1FGL J0005.7+3815", "Check source name");
    test_value(sky5.instruments(), "LAT,CTA", "Check instruments");
    test_value(sky5.ids(), "000001", "Check observation identifier");
    test_value(sky5.ts(), 100.1, 1.0e-6, "Check TS value");
    test_assert(sky5.tscalc(), "Check whether TS computation is requested");
    test_assert(sky5.has_scales(), "Check whether sky model has scales");
    test_assert(sky5.has_ts(), "Check whether sky model has TS value");
    test_value(sky5.scale("LAT").value(), 1.1, 1.0e-6, "Check LAT scale factor");
    test_value(sky5.scale("CTA").value(), 0.5, 1.0e-6, "Check CTA scale factor");
    test_value(sky5.type(), "PointSource", "Check source type");
    test_assert(sky5.spatial()  != NULL, "Check whether sky model has spatial component");
    test_assert(sky5.spectral() != NULL, "Check whether sky model has spectral component");
    test_assert(sky5.temporal() != NULL, "Check whether sky model has temporal component");

    // Test read method
    sky5.clear();
    sky5.read(*element1);
    test_value(sky5.size(), 6, "Check number of model parameters");
    test_value(sky5.name(), "1FGL J0005.7+3815", "Check source name");
    test_value(sky5.instruments(), "LAT,CTA", "Check instruments");
    test_value(sky5.ids(), "000001", "Check observation identifier");
    test_value(sky5.ts(), 100.1, 1.0e-6, "Check TS value");
    test_assert(sky5.tscalc(), "Check whether TS computation is requested");
    test_assert(sky5.has_scales(), "Check whether sky model has scales");
    test_assert(sky5.has_ts(), "Check whether sky model has TS value");
    test_value(sky5.scale("LAT").value(), 1.1, 1.0e-6, "Check LAT scale factor");
    test_value(sky5.scale("CTA").value(), 0.5, 1.0e-6, "Check CTA scale factor");
    test_value(sky5.type(), "PointSource", "Check source type");
    test_assert(sky5.spatial()  != NULL, "Check whether sky model has spatial component");
    test_assert(sky5.spectral() != NULL, "Check whether sky model has spectral component");
    test_assert(sky5.temporal() != NULL, "Check whether sky model has temporal component");

    // Test model writing and back reading
    GModels models;
    models.append(sky5);
    models.save("test_sky_model.xml");
    models.clear();
    models.load("test_sky_model.xml");
    sky5.clear();
    sky5 = *(static_cast<GModelSky*>(models[0]));
    test_value(sky5.size(), 6, "Check number of model parameters");
    test_value(sky5.name(), "1FGL J0005.7+3815", "Check source name");
    test_value(sky5.instruments(), "LAT,CTA", "Check instruments");
    test_value(sky5.ids(), "000001", "Check observation identifier");
    test_value(sky5.ts(), 100.1, 1.0e-6, "Check TS value");
    test_assert(sky5.tscalc(), "Check whether TS computation is requested");
    test_assert(sky5.has_scales(), "Check whether sky model has scales");
    test_assert(sky5.has_ts(), "Check whether sky model has TS value");
    test_value(sky5.scale("LAT").value(), 1.1, 1.0e-6, "Check LAT scale factor");
    test_value(sky5.scale("CTA").value(), 0.5, 1.0e-6, "Check CTA scale factor");
    test_value(sky5.type(), "PointSource", "Check source type");
    test_assert(sky5.spatial()  != NULL, "Check whether sky model has spatial component");
    test_assert(sky5.spectral() != NULL, "Check whether sky model has spectral component");
    test_assert(sky5.temporal() != NULL, "Check whether sky model has temporal component");

    // Test value method
    GSkyDir dir;
    dir.radec_deg(83.6331, 22.0145);
    GEnergy energy(100.0, "MeV");
    GTime   time(0.0);
    test_value(sky5.value(GPhoton(dir, energy, time)), 1.73e-07);

    // Test gradient method
    GVector vector = sky5.gradients(GPhoton(dir, energy, time));
    test_value(vector[0], 0.0);
    test_value(vector[1], 0.0);
    test_value(vector[2], 1e-07);
    test_value(vector[3], 0.0);
    test_value(vector[4], 0.0);
    test_value(vector[5], 0.0);

    // Get test sky model
    GXml         xml2(m_xml_file);
    GXmlElement* element2 = xml2.element(0)->element(0);
    GModelSky    sky6(*element2);

    // Test has_par methods
    test_assert(sky6.spatial()->has_par("RA"),
                "Expect \"RA\" parameter in spatial component");
    test_assert(!sky6.spatial()->has_par("RAX"),
                "Do not expect \"RAX\" parameter in spatial component");
    test_assert(sky6.spectral()->has_par("Prefactor"),
                "Expect \"Prefactor\" parameter in spectral component");
    test_assert(!sky6.spectral()->has_par("Prefactors"),
                "Do not expect \"Prefactors\" parameter in spectral component");
    test_assert(sky6.temporal()->has_par("Normalization"),
                "Expect \"Normalization\" parameter in temporal component");
    test_assert(!sky6.temporal()->has_par("Normalizations"),
                "Do not expect \"Normalizations\" parameter in temporal component");

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Test GModelSpatialPointSource class
 ***************************************************************************/
void TestGModel::test_point_source(void)
{
    // Test void constructor
    test_try("Test void constructor");
    try {
        GModelSpatialPointSource model;
        test_assert(model.type() == "PointSource",
                                    "Model type \"PointSource\" expected.");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test sky direction constructor
    test_try("Test sky direction constructor");
    try {
        GSkyDir dir;
        dir.radec_deg(83.6331, +22.0145);
        GModelSpatialPointSource model(dir);
        test_value(model.ra(), 83.6331);
        test_value(model.dec(), +22.0145);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test value constructor
    test_try("Test value constructor");
    try {
        GModelSpatialPointSource model(83.6331, +22.0145);
        test_value(model.ra(), 83.6331);
        test_value(model.dec(), +22.0145);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }
    
    // Test XML constructor and value
    test_try("Test XML constructor, value and gradients");
    try {
        // Test XML constructor
        GXml                     xml(m_xml_model_point_plaw);
        GXmlElement*             element = xml.element(0)->element(0)->element("spatialModel", 0);
        GModelSpatialPointSource model(*element);
        test_value(model.size(), 2);
        test_assert(model.type() == "PointSource", "Expected \"PointSource\"");
        test_value(model.ra(), 83.6331);
        test_value(model.dec(), +22.0145);

        // Test ra method
        model.ra(3.9);
        test_value(model.ra(), 3.9);

        // Test dec method
        model.dec(3.9);
        test_value(model.dec(), 3.9);

        // Test operator access
        const char* strarray[] = {"RA", "DEC"};
        for (int i = 0; i < 2; ++i) {
            std::string keyname(strarray[i]);
            model[keyname].value(2.1);
            model[keyname].error(1.9);
            model[keyname].gradient(0.8);
            test_value(model[keyname].value(), 2.1);
            test_value(model[keyname].error(), 1.9);
            test_value(model[keyname].gradient(), 0.8);
        }

        // Success if we reached this point
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Test GModelSpatialDiffuseConst class
 ***************************************************************************/
void TestGModel::test_diffuse_const(void)
{
    // Test void constructor
    test_try("Test void constructor");
    try {
        GModelSpatialDiffuseConst model;
        test_assert(model.type() == "DiffuseIsotropic",
                                    "Model type \"DiffuseIsotropic\" expected.");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test value constructor
    test_try("Test value constructor");
    try {
        GModelSpatialDiffuseConst model(3.0);
        test_value(model.value(), 3.0);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }
    
    // Test XML constructor and value
    test_try("Test XML constructor, value and gradients");
    try {
        // Test XML constructor
        GXml                      xml(m_xml_model_diffuse_const);
        GXmlElement*              element = xml.element(0)->element(0)->element("spatialModel", 0);
        GModelSpatialDiffuseConst model(*element);
        test_value(model.size(), 1);
        test_assert(model.type() == "DiffuseIsotropic",
                    "Expected \"DiffuseIsotropic\"");
        test_value(model.value(), 1.0);

        // Test value method
        model.value(3.9);
        test_value(model.value(), 3.9);

        // Test operator access
        test_value(model["Value"].value(), 3.9);
        test_value(model["Value"].error(), 0.0);
        test_value(model["Value"].gradient(), 0.0);
        model["Value"].value(2.1);
        model["Value"].error(1.9);
        model["Value"].gradient(0.8);
        test_value(model["Value"].value(), 2.1);
        test_value(model["Value"].error(), 1.9);
        test_value(model["Value"].gradient(), 0.8);

        // Success if we reached this point
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Test GModelSpatialDiffuseCube class
 ***************************************************************************/
void TestGModel::test_diffuse_cube(void)
{
    // Test void constructor
    test_try("Test void constructor");
    try {
        GModelSpatialDiffuseCube model;
        test_assert(model.type() == "DiffuseMapCube",
                                    "Model type \"DiffuseMapCube\" expected.");
        test_assert(model.filename().url() == "", "Model filename \"\" expected.");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test filename value constructor
    test_try("Test filename value constructor");
    try {
        GModelSpatialDiffuseCube model(m_cube_file, 3.0);
        test_value(model.value(), 3.0);
        test_assert(model.filename().url() == m_cube_file,
                    "Expected \""+m_cube_file+"\"");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test skymap value constructor
    test_try("Test skymap value constructor");
    try {
        GSkyMap map("GAL", 16, "RING", 10);
        GEnergies energies;
        for (int i = 0; i < 10; ++i) {
            energies.append(GEnergy(double(i+1.0), "MeV"));
        }
        GModelSpatialDiffuseCube model(map, energies, 3.0);
        test_value(model.value(), 3.0);
        test_assert(model.filename().url() == "", "Expected \"\"");
        //test_assert(model.cube() == map, "Map cube is not the expected one");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }
    
    // Test XML constructor and attribute methods
    test_try("Test XML constructor, value and gradients");
    try {
        // Test XML constructor
        GXml                     xml(m_xml_model_diffuse_cube);
        GXmlElement*             element = xml.element(0)->element(0)->element("spatialModel", 0);
        GModelSpatialDiffuseCube model(*element);
        test_value(model.size(), 1);
        test_assert(model.type() == "DiffuseMapCube",
                    "Expected \"DiffuseMapCube\"");
        test_value(model.value(), 1.0);
        test_assert(model.filename().url() == m_cube_file,
                    "Model filename \""+m_cube_file+"\" expected.");

        // Test value method
        model.value(3.9);
        test_value(model.value(), 3.9);

        // Test filename method
        model.filename("Help me!");
        test_assert(model.filename().url() == "Help me!",
                    "Model filename \"Help me!\" expected.");

        // Test cube method
        model.cube(GSkyMap("GAL", 16, "RING", 10));
        test_value(model.cube().npix(), 3072);

        // Test operator access
        test_value(model["Normalization"].value(), 3.9);
        test_value(model["Normalization"].error(), 0.0);
        test_value(model["Normalization"].gradient(), 0.0);
        model["Normalization"].value(2.1);
        model["Normalization"].error(1.9);
        model["Normalization"].gradient(0.8);
        test_value(model["Normalization"].value(), 2.1);
        test_value(model["Normalization"].error(), 1.9);
        test_value(model["Normalization"].gradient(), 0.8);

        // Success if we reached this point
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Test GModelSpatialDiffuseMap class
 ***************************************************************************/
void TestGModel::test_diffuse_map(void)
{
    // Test void constructor
    GModelSpatialDiffuseMap model1;
    test_assert(model1.type() == "DiffuseMap",
                "Check that model is of type \"DiffuseMap\".");
    test_assert(model1.filename().url() == "",
                "Check that model has empty filename.");

    // Test filename value constructor
    GModelSpatialDiffuseMap model2(m_map_file, 3.0);
    test_value(model2.value(), 3.0,
               "Check the value after loading from XML file.");
    test_assert(model2.filename().url() == m_map_file,
                "Check the filename after loading from XML file.");

    // Test skymap value constructor
    GSkyMap map("GAL", 16, "RING", 10);
    GModelSpatialDiffuseMap model3(map, 3.0);
    test_value(model3.value(), 3.0,
               "Check the value after creating model from sky map.");
    test_assert(model3.filename().url() == "",
                "Check that filename is empty after creating model from sky map.");
    
    // Test XML constructor and attribute methods
    GXml         xml(m_xml_model_diffuse_map);
    GXmlElement* elementp = xml.element(0)->element(0)->element("spatialModel", 0);
    GModelSpatialDiffuseMap model(*elementp);
    test_value(model.size(), 1,
               "Check that model build from an XML element has one parameter.");
    test_assert(model.type() == "DiffuseMap",
                "Check that model build from an XML element is of type "
                "\"DiffuseMap\".");
    test_value(model.value(), 1.0,
               "Check that normalisation of model build from an XML element is "
               "unity.");
    test_assert(model.filename().url() == m_map_file,
                "Check that filename of model build from an XML elements is "
                "\""+m_map_file+"\".");

    // Test value method
    model.value(3.9);
    test_value(model.value(), 3.9);

    // Test load method
    model.load(m_map_file);
    test_assert(model.filename().url() == m_map_file,
                "Check that filename of model loaded from a file is "
                "\""+m_map_file+"\".");

    // Test map method
    model.map(GSkyMap("GAL", 16, "RING", 10));
    test_value(model.map().npix(), 3072);

    // Test operator access
    test_value(model["Normalization"].value(), 3.9);
    test_value(model["Normalization"].error(), 0.0);
    test_value(model["Normalization"].gradient(), 0.0);
    model["Normalization"].value(2.1);
    model["Normalization"].error(1.9);
    model["Normalization"].gradient(0.8);
    test_value(model["Normalization"].value(), 2.1);
    test_value(model["Normalization"].error(), 1.9);
    test_value(model["Normalization"].gradient(), 0.8);

    // Initialize photon direction for testing
    GSkyDir dir;
    GEnergy energy;
    GTime   time;
    dir.radec_deg(201.68111, -42.68742);
    GPhoton photon(dir, energy, time);

    // Test normalized map
    GModelSpatialDiffuseMap map_norm(m_map_file, 3.0, true);
    test_assert(map_norm.normalize(),
                "Check that model is normalized.");
    test_value(map_norm.eval(photon), 13069.7741687, 0.2,
               "Check that model has the correct skymap intensity.");
    GXmlElement element;
    map_norm.write(element);
    map_norm.read(element);
    test_assert(map_norm.normalize(),
                "Check that model is normalized.");
    test_value(map_norm.eval(photon), 13069.7741687, 0.2,
               "Check that model has the correct skymap intensity.");

    // Test non-normalized map
    GModelSpatialDiffuseMap map_nonnorm(m_map_file, 3.0, false);
    test_assert(!map_nonnorm.normalize(),
                "Check that model is not normalized.");
    test_value(map_nonnorm.eval(photon), 13069.6002755, 0.2,
               "Check that model has the correct skymap intensity.");
    GXmlElement element2;
    map_nonnorm.write(element2);
    map_nonnorm.read(element2);
    test_assert(!map_nonnorm.normalize(),
                "Check that model is not normalized.");

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Test GModelSpatialRadialDisk class
 ***************************************************************************/
void TestGModel::test_radial_disk(void)
{
    // Test void constructor
    test_try("Test void constructor");
    try {
        GModelSpatialRadialDisk model;
        test_assert(model.type() == "RadialDisk",
                                    "Model type \"RadialDisk\" expected.");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test value constructor
    test_try("Test value constructor");
    try {
        GSkyDir dir;
        dir.radec_deg(83.6331, +22.0145);
        GModelSpatialRadialDisk model(dir, 3.0);
        test_value(model.ra(), 83.6331);
        test_value(model.dec(), 22.0145);
        test_value(model.radius(), 3.0);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }
    
    // Test XML constructor and attribute methods
    test_try("Test XML constructor, value and gradients");
    try {
        // Test XML constructor
        GXml                    xml(m_xml_model_radial_disk);
        GXmlElement*            element = xml.element(0)->element(0)->element("spatialModel", 0);
        GModelSpatialRadialDisk model(*element);
        test_value(model.size(), 3);
        test_assert(model.type() == "RadialDisk", "Expected \"RadialDisk\"");
        test_value(model.ra(), 83.6331);
        test_value(model.dec(), 22.0145);
        test_value(model.radius(), 0.45);

        // Test ra method
        model.ra(100.0);
        test_value(model.ra(), 100.0);

        // Test dec method
        model.dec(10.0);
        test_value(model.dec(), 10.0);

        // Test dir method
        GSkyDir dir;
        dir.radec_deg(83.6331, +22.0145);
        model.dir(dir);
        test_assert(model.dir() == dir, "Test sky direction");

        // Test radius method
        model.radius(3.9);
        test_value(model.radius(), 3.9);

        // Test operator access
        const char* strarray[] = {"RA", "DEC", "Radius"};
        for (int i = 0; i < 3; ++i) {
            std::string keyname(strarray[i]);
            model[keyname].value(2.1);
            model[keyname].error(1.9);
            model[keyname].gradient(0.8);
            test_value(model[keyname].value(), 2.1);
            test_value(model[keyname].error(), 1.9);
            test_value(model[keyname].gradient(), 0.8);
        }

        // Success if we reached this point
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Test GModelSpatialRadialGauss class
 ***************************************************************************/
void TestGModel::test_radial_gauss(void)
{
    // Test void constructor
    test_try("Test void constructor");
    try {
        GModelSpatialRadialGauss model;
        test_assert(model.type() == "RadialGaussian",
                                    "Model type \"RadialGaussian\" expected.");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test value constructor
    test_try("Test value constructor");
    try {
        GSkyDir dir;
        dir.radec_deg(83.6331, +22.0145);
        GModelSpatialRadialGauss model(dir, 3.0);
        test_value(model.ra(), 83.6331);
        test_value(model.dec(), 22.0145);
        test_value(model.sigma(), 3.0);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }
    
    // Test XML constructor and attribute methods
    test_try("Test XML constructor, value and gradients");
    try {
        // Test XML constructor
        GXml                     xml(m_xml_model_radial_gauss);
        GXmlElement*             element = xml.element(0)->element(0)->element("spatialModel", 0);
        GModelSpatialRadialGauss model(*element);
        test_value(model.size(), 3);
        test_assert(model.type() == "RadialGaussian",
                    "Expected \"RadialGaussian\"");
        test_value(model.ra(), 83.6331);
        test_value(model.dec(), 22.0145);
        test_value(model.sigma(), 0.20);

        // Test ra method
        model.ra(100.0);
        test_value(model.ra(), 100.0);

        // Test dec method
        model.dec(10.0);
        test_value(model.dec(), 10.0);

        // Test dir method
        GSkyDir dir;
        dir.radec_deg(83.6331, +22.0145);
        model.dir(dir);
        test_assert(model.dir() == dir, "Test sky direction");

        // Test sigma method
        model.sigma(3.9);
        test_value(model.sigma(), 3.9);

        // Test operator access
        const char* strarray[] = {"RA", "DEC", "Sigma"};
        for (int i = 0; i < 3; ++i) {
            std::string keyname(strarray[i]);
            model[keyname].value(2.1);
            model[keyname].error(1.9);
            model[keyname].gradient(0.8);
            test_value(model[keyname].value(), 2.1);
            test_value(model[keyname].error(), 1.9);
            test_value(model[keyname].gradient(), 0.8);
        }

        // Success if we reached this point
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Test GModelSpatialRadialShell class
 ***************************************************************************/
void TestGModel::test_radial_shell(void)
{
    // Test void constructor
    test_try("Test void constructor");
    try {
        GModelSpatialRadialShell model;
        test_assert(model.type() == "RadialShell",
                                    "Model type \"RadialShell\" expected.");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test value constructor
    test_try("Test value constructor");
    try {
        GSkyDir dir;
        dir.radec_deg(83.6331, +22.0145);
        GModelSpatialRadialShell model(dir, 3.0, 1.0);
        test_value(model.ra(), 83.6331);
        test_value(model.dec(), 22.0145);
        test_value(model.radius(), 3.0);
        test_value(model.width(), 1.0);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }
    
    // Test XML constructor and attribute methods
    test_try("Test XML constructor, value and gradients");
    try {
        // Test XML constructor
        GXml                     xml(m_xml_model_radial_shell);
        GXmlElement*             element = xml.element(0)->element(0)->element("spatialModel", 0);
        GModelSpatialRadialShell model(*element);
        test_value(model.size(), 4);
        test_assert(model.type() == "RadialShell", "Expected \"RadialShell\"");
        test_value(model.ra(), 83.6331);
        test_value(model.dec(), 22.0145);
        test_value(model.radius(), 0.30);
        test_value(model.width(), 0.10);

        // Test ra method
        model.ra(100.0);
        test_value(model.ra(), 100.0);

        // Test dec method
        model.dec(10.0);
        test_value(model.dec(), 10.0);

        // Test dir method
        GSkyDir dir;
        dir.radec_deg(83.6331, +22.0145);
        model.dir(dir);
        test_assert(model.dir() == dir, "Test sky direction");

        // Test radius method
        model.radius(3.9);
        test_value(model.radius(), 3.9);

        // Test width method
        model.width(3.9);
        test_value(model.width(), 3.9);

        // Test operator access
        const char* strarray[] = {"RA", "DEC", "Radius", "Width"};
        for (int i = 0; i < 4; ++i) {
            std::string keyname(strarray[i]);
            model[keyname].value(2.1);
            model[keyname].error(1.9);
            model[keyname].gradient(0.8);
            test_value(model[keyname].value(), 2.1);
            test_value(model[keyname].error(), 1.9);
            test_value(model[keyname].gradient(), 0.8);
        }

        // Success if we reached this point
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Test GModelSpatialEllipticalDisk class
 ***************************************************************************/
void TestGModel::test_elliptical_disk(void)
{
    // Test void constructor
    test_try("Test void constructor");
    try {
        GModelSpatialEllipticalDisk model;
        test_assert(model.type() == "EllipticalDisk",
                                    "Model type \"EllipticalDisk\" expected.");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test value constructor
    test_try("Test value constructor");
    try {
        GSkyDir dir;
        dir.radec_deg(83.6331, +22.0145);
        GModelSpatialEllipticalDisk model(dir, 3.0, 2.0, 45.0);
        test_value(model.ra(), 83.6331);
        test_value(model.dec(), 22.0145);
        test_value(model.posangle(), 45.0);
        test_value(model.semimajor(), 3.0);
        test_value(model.semiminor(), 2.0);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }
    
    // Test XML constructor and attribute methods
    test_try("Test XML constructor, value and gradients");
    try {
        // Test XML constructor
        GXml                        xml(m_xml_model_elliptical_disk);
        GXmlElement*                element = xml.element(0)->element(0)->element("spatialModel", 0);
        GModelSpatialEllipticalDisk model(*element);
        test_value(model.size(), 5);
        test_assert(model.type() == "EllipticalDisk", "Expected \"EllipticalDisk\"");
        test_value(model.ra(), 83.6331);
        test_value(model.dec(), 22.0145);
        test_value(model.posangle(), 45.0);
        test_value(model.semimajor(), 2.0);
        test_value(model.semiminor(), 0.5);

        // Test ra method
        model.ra(100.0);
        test_value(model.ra(), 100.0);

        // Test dec method
        model.dec(10.0);
        test_value(model.dec(), 10.0);

        // Test dir method
        GSkyDir dir;
        dir.radec_deg(83.6331, +22.0145);
        model.dir(dir);
        test_assert(model.dir() == dir, "Test sky direction");

        // Test posangle method
        model.posangle(3.9);
        test_value(model.posangle(), 3.9);

        // Test semimajor method
        model.semimajor(3.9);
        test_value(model.semimajor(), 3.9);

        // Test semiminor method
        model.semiminor(3.9);
        test_value(model.semiminor(), 3.9);

        // Test operator access
        const char* strarray[] = {"RA", "DEC", "PA", "MinorRadius", "MajorRadius"};
        for (int i = 0; i < 5; ++i) {
            std::string keyname(strarray[i]);
            model[keyname].value(2.1);
            model[keyname].error(1.9);
            model[keyname].gradient(0.8);
            test_value(model[keyname].value(), 2.1);
            test_value(model[keyname].error(), 1.9);
            test_value(model[keyname].gradient(), 0.8);
        }

        // Success if we reached this point
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Exit test
    return;
}

/***********************************************************************//**
 * @brief Test GModelSpatialEllipticalGauss class
 ***************************************************************************/
void TestGModel::test_elliptical_gauss(void)
{
    // Test void constructor
    test_try("Test void constructor");
    try {
        GModelSpatialEllipticalGauss model;
        test_assert(model.type() == "EllipticalGaussian",
                                    "Model type \"EllipticalGaussian\" expected.");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test value constructor
    test_try("Test value constructor");
    try {
        GSkyDir dir;
        dir.radec_deg(83.6331, +22.0145);
        GModelSpatialEllipticalGauss model(dir, 3.0, 2.0, 45.0);
        test_value(model.ra(), 83.6331);
        test_value(model.dec(), 22.0145);
        test_value(model.posangle(), 45.0);
        test_value(model.semimajor(), 3.0);
        test_value(model.semiminor(), 2.0);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test XML constructor and attribute methods
    test_try("Test XML constructor, value and gradients");
    try {
        // Test XML constructor
        GXml                         xml(m_xml_model_elliptical_gauss);
        GXmlElement*                 element = xml.element(0)->element(0)->element("spatialModel", 0);
        GModelSpatialEllipticalGauss model(*element);
        test_value(model.size(), 5);
        test_assert(model.type() == "EllipticalGaussian",
                    "Expected \"EllipticalGaussian\"");
        test_value(model.ra(), 83.6331);
        test_value(model.dec(), 22.0145);
        test_value(model.posangle(), 45.0);
        test_value(model.semimajor(), 0.3);
        test_value(model.semiminor(), 0.1);

        // Test ra method
        model.ra(100.0);
        test_value(model.ra(), 100.0);

        // Test dec method
        model.dec(10.0);
        test_value(model.dec(), 10.0);

        // Test dir method
        GSkyDir dir;
        dir.radec_deg(83.6331, +22.0145);
        model.dir(dir);
        test_assert(model.dir() == dir, "Test sky direction");

        // Test posangle method
        model.posangle(3.9);
        test_value(model.posangle(), 3.9);

        // Test semimajor method
        model.semimajor(3.9);
        test_value(model.semimajor(), 3.9);

        // Test semiminor method
        model.semiminor(3.9);
        test_value(model.semiminor(), 3.9);

        // Test operator access
        const char* strarray[] = {"RA", "DEC", "PA", "MinorRadius", "MajorRadius"};
        for (int i = 0; i < 5; ++i) {
            std::string keyname(strarray[i]);
            model[keyname].value(2.1);
            model[keyname].error(1.9);
            model[keyname].gradient(0.8);
            test_value(model[keyname].value(), 2.1);
            test_value(model[keyname].error(), 1.9);
            test_value(model[keyname].gradient(), 0.8);
        }

        // Success if we reached this point
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Exit test
    return;
}



/***********************************************************************//**
 * @brief Test GModelSpectralConst class
 ***************************************************************************/
void TestGModel::test_const(void)
{
    // Test void constructor
    GModelSpectralConst model1;
    test_value(model1.type(), "Constant", "Check void model type");

    // Test value constructor
    GModelSpectralConst model2(3.0);
    test_value(model2.value(), 3.0);
    
    // Test XML constructor
    GXml                xml(m_xml_model_point_const);
    GXmlElement*        element = xml.element(0)->element(0)->element("spectrum", 0);
    GModelSpectralConst model3(*element);
    test_value(model3.size(), 1);
    test_value(model3.type(), "Constant", "Check model type");
    test_value(model3.value(), 5.7e-16);

    // Test value method
    model3.value(3.9e-16);
    test_value(model3.value(), 3.9e-16);

    // Test operator access
    const char* strarray[] = {"Normalization"};
    for (int i = 0; i < 1; ++i) {
        std::string keyname(strarray[i]);
        model3[keyname].remove_range(); // To allow setting of any value
        model3[keyname].value(2.1);
        model3[keyname].error(1.9);
        model3[keyname].gradient(0.8);
        test_value(model3[keyname].value(), 2.1);
        test_value(model3[keyname].error(), 1.9);
        test_value(model3[keyname].gradient(), 0.8);
    }

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Test GModelSpectralGauss class
 ***************************************************************************/
void TestGModel::test_gauss(void)
{
    // Test void constructor
    GModelSpectralGauss model1;
    test_value(model1.type(), "Gaussian", "Check void model type");

    // Test value constructor
    GModelSpectralGauss model2(1.0, GEnergy(42.0, "MeV"), GEnergy(0.5, "MeV"));
    test_value(model2.norm(), 1.0);
    test_value(model2.mean().MeV(), 42.0);
    test_value(model2.sigma().MeV(), 0.5);

    // Test XML constructor
    GXml                xml(m_xml_model_point_gauss);
    GXmlElement*        element = xml.element(0)->element(0)->element("spectrum", 0);
    GModelSpectralGauss model3(*element);
    test_value(model3.size(), 3);
    test_value(model3.type(), "Gaussian", "Check model type");
    test_value(model3.norm(), 1.0e-10);
    test_value(model3.mean().MeV(), 5.0e6);
    test_value(model3.sigma().MeV(), 1.0e6);

    // Test norm method
    model3.norm(2.3e-10);
    test_value(model3.norm(), 2.3e-10);

    // Test mean method
    model3.mean(GEnergy(1.5678, "TeV"));
    test_value(model3.mean().TeV(), 1.5678);

    // Test sigma method
    model3.sigma(GEnergy(0.1234, "TeV"));
    test_value(model3.sigma().TeV(), 0.1234);

    // Test operator access
    const char* strarray[] = {"Normalization", "Mean", "Sigma"};
    for (int i = 0; i < 3; ++i) {
        std::string keyname(strarray[i]);
        model3[keyname].remove_range(); // To allow setting of any value
        model3[keyname].value(2.1);
        model3[keyname].error(1.9);
        model3[keyname].gradient(0.8);
        test_value(model3[keyname].value(), 2.1);
        test_value(model3[keyname].error(), 1.9);
        test_value(model3[keyname].gradient(), 0.8);
    }

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Test GModelSpectralPlaw class
 ***************************************************************************/
void TestGModel::test_plaw(void)
{
    // Test void constructor
    GModelSpectralPlaw model1;
    test_value(model1.type(), "PowerLaw", "Check type of void model");

    // Test value constructor
    GModelSpectralPlaw model2(2.0, -2.1, GEnergy(100.0, "MeV"));
    test_value(model2.prefactor(), 2.0);
    test_value(model2.index(), -2.1);
    test_value(model2.pivot().MeV(), 100.0);
    
    // Test XML constructor
    GXml               xml(m_xml_model_point_plaw);
    GXmlElement*       element = xml.element(0)->element(0)->element("spectrum", 0);
    GModelSpectralPlaw model3(*element);
    test_value(model3.size(), 3);
    test_value(model3.type(), "PowerLaw", "Check model type");
    test_value(model3.prefactor(), 5.7e-16);
    test_value(model3.index(), -2.48);
    test_value(model3.pivot().TeV(), 0.3);

    // Test prefactor method
    model3["Prefactor"].remove_range(); // To allow setting of any value
    model3.prefactor(3.9);
    test_value(model3.prefactor(), 3.9);

    // Test index method
    model3["Index"].remove_range(); // To allow setting of any value
    model3.index(2.1);
    test_value(model3.index(), 2.1);

    // Test pivot method
    model3["PivotEnergy"].remove_range(); // To allow setting of any value
    model3.pivot(GEnergy(10.0, "MeV"));
    test_value(model3.pivot().MeV(), 10.0);

    // Test operator access
    const char* strarray[] = {"Prefactor", "Index", "PivotEnergy"};
    for (int i = 0; i < 3; ++i) {
        std::string keyname(strarray[i]);
        model3[keyname].value(2.1);
        model3[keyname].error(1.9);
        model3[keyname].gradient(0.8);
        test_value(model3[keyname].value(), 2.1);
        test_value(model3[keyname].error(), 1.9);
        test_value(model3[keyname].gradient(), 0.8);
    }

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Test GModelSpectralPlawPhotonFlux class
 ***************************************************************************/
void TestGModel::test_plaw_phflux(void)
{
    // Test void constructor
    GModelSpectralPlawPhotonFlux model1;
    test_value(model1.type(), "PowerLaw", "Check type of void model");

    // Test value constructor
    GModelSpectralPlawPhotonFlux model2(2.0, -2.1, GEnergy(10.0, "MeV"),
                                                   GEnergy(100.0, "MeV"));
    test_value(model2.flux(), 2.0);
    test_value(model2.index(), -2.1);
    test_value(model2.emin().MeV(), 10.0);
    test_value(model2.emax().MeV(), 100.0);

    // Test XML constructor and value
    GXml         xml(m_xml_model_point_plaw_phflux);
    GXmlElement* element = xml.element(0)->element(0)->element("spectrum", 0);
    GModelSpectralPlawPhotonFlux model3(*element);
    test_value(model3.size(), 4);
    test_value(model3.type(), "PowerLaw", "Check model type");
    test_value(model3.flux(), 1.0e-7);
    test_value(model3.index(), -2.0);
    test_value(model3.emin().MeV(), 100.0);
    test_value(model3.emax().MeV(), 500000.0);

    // Test integral method
    model3.flux(2.1e-7);
    test_value(model3.flux(), 2.1e-7);

    // Test index method
    model3.index(-2.3);
    test_value(model3.index(), -2.3);

    // Test emin method
    model3.emin(GEnergy(10.0, "MeV"));
    test_value(model3.emin().MeV(), 10.0);

    // Test emax method
    model3.emax(GEnergy(10.0, "MeV"));
    test_value(model3.emax().MeV(), 10.0);

    // Test operator access
    const char* strarray[] = {"PhotonFlux", "Index", "LowerLimit", "UpperLimit"};
    for (int i = 0; i < 4; ++i) {
        std::string keyname(strarray[i]);
        model3[keyname].remove_range(); // To allow setting of any value
        model3[keyname].value(2.1);
        model3[keyname].error(1.9);
        model3[keyname].gradient(0.8);
        test_value(model3[keyname].value(), 2.1);
        test_value(model3[keyname].error(), 1.9);
        test_value(model3[keyname].gradient(), 0.8);
    }

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Test GModelSpectralPlawEnergyFlux class
 ***************************************************************************/
void TestGModel::test_plaw_eflux(void)
{
    // Test void constructor
    GModelSpectralPlawEnergyFlux model1;
    test_value(model1.type(), "PowerLaw", "Check type of void model");

    // Test value constructor
    GModelSpectralPlawEnergyFlux model2(2.0, -2.1, GEnergy(10.0, "MeV"),
                                                   GEnergy(100.0, "MeV"));
    test_value(model2.eflux(), 2.0);
    test_value(model2.index(), -2.1);
    test_value(model2.emin().MeV(), 10.0);
    test_value(model2.emax().MeV(), 100.0);
    
    // Test XML constructor and value
    GXml         xml(m_xml_model_point_plaw_eflux);
    GXmlElement* element = xml.element(0)->element(0)->element("spectrum", 0);
    GModelSpectralPlawEnergyFlux model3(*element);
    test_value(model3.size(), 4);
    test_value(model3.type(), "PowerLaw", "Check model type");
    test_value(model3.eflux(), 1.0e-7);
    test_value(model3.index(), -2.0);
    test_value(model3.emin().MeV(), 100.0);
    test_value(model3.emax().MeV(), 500000.0);

    // Test integral method
    model3.eflux(2.1e-7);
    test_value(model3.eflux(), 2.1e-7);

    // Test index method
    model3.index(-2.3);
    test_value(model3.index(), -2.3);

    // Test emin method
    model3.emin(GEnergy(10.0, "MeV"));
    test_value(model3.emin().MeV(), 10.0);

    // Test emax method
    model3.emax(GEnergy(10.0, "MeV"));
    test_value(model3.emax().MeV(), 10.0);

    // Test operator access
    const char* strarray[] = {"EnergyFlux", "Index", "LowerLimit", "UpperLimit"};
    for (int i = 0; i < 4; ++i) {
        std::string keyname(strarray[i]);
        model3[keyname].remove_range(); // To allow setting of any value
        model3[keyname].value(2.1);
        model3[keyname].error(1.9);
        model3[keyname].gradient(0.8);
        test_value(model3[keyname].value(), 2.1);
        test_value(model3[keyname].error(), 1.9);
        test_value(model3[keyname].gradient(), 0.8);
    }

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Test GModelSpectralExpPlaw class
 ***************************************************************************/
void TestGModel::test_eplaw(void)
{
    // Test void constructor
    GModelSpectralExpPlaw model1;
    test_value(model1.type(), "ExponentialCutoffPowerLaw", "Check type of void model");

    // Test value constructor
    GModelSpectralExpPlaw model2(2.0, -2.1, GEnergy(100.0, "MeV"), GEnergy(1.0, "GeV"));
    test_value(model2.prefactor(), 2.0);
    test_value(model2.index(), -2.1);
    test_value(model2.pivot().MeV(), 100.0);
    test_value(model2.cutoff().GeV(), 1.0);
    
    // Test XML constructor
    GXml                  xml(m_xml_model_point_eplaw);
    GXmlElement*          element = xml.element(0)->element(0)->element("spectrum", 0);
    GModelSpectralExpPlaw model3(*element);
    test_value(model3.size(), 4);
    test_value(model3.type(), "ExponentialCutoffPowerLaw", "Check model type");
    test_value(model3.prefactor(), 5.7e-16);
    test_value(model3.index(), -2.48);
    test_value(model3.pivot().TeV(), 0.3);
    test_value(model3.cutoff().TeV(), 1.0);

    // Test prefactor method
    model3.prefactor(2.3e-16);
    test_value(model3.prefactor(), 2.3e-16);

    // Test index method
    model3.index(-2.6);
    test_value(model3.index(), -2.6);

    // Test pivot method
    model3.pivot(GEnergy(0.5, "TeV"));
    test_value(model3.pivot().TeV(), 0.5);

    // Test cutoff method
    model3.cutoff(GEnergy(10.0, "TeV"));
    test_value(model3.cutoff().TeV(), 10.0);

    // Test operator access
    const char* strarray[] = {"Prefactor", "Index", "PivotEnergy", "CutoffEnergy"};
    for (int i = 0; i < 4; ++i) {
        std::string keyname(strarray[i]);
        model3[keyname].remove_range(); // To allow setting of any value
        model3[keyname].value(2.1);
        model3[keyname].error(1.9);
        model3[keyname].gradient(0.8);
        test_value(model3[keyname].value(), 2.1);
        test_value(model3[keyname].error(), 1.9);
        test_value(model3[keyname].gradient(), 0.8);
    }

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Test GModelSpectralExpInvPlaw class
 ***************************************************************************/
void TestGModel::test_einvplaw(void)
{
    // Test void constructor
    GModelSpectralExpInvPlaw model1;
    test_value(model1.type(), "ExponentialCutoffPowerLaw", "Check type of void model");

    // Test value constructor
    GModelSpectralExpInvPlaw model2(2.0, -2.1, GEnergy(100.0, "MeV"), 1.0e-3);
    test_value(model2.prefactor(), 2.0);
    test_value(model2.index(), -2.1);
    test_value(model2.pivot().MeV(), 100.0);
    test_value(model2.inverse_cutoff(), 1.0e-3);
    test_value(model2.cutoff().GeV(), 1.0);

    // Test alternative value constructor
    GModelSpectralExpInvPlaw model3(2.0, -2.1, GEnergy(100.0, "MeV"),
                                               GEnergy(1.0, "GeV"));
    test_value(model3.prefactor(), 2.0);
    test_value(model3.index(), -2.1);
    test_value(model3.pivot().MeV(), 100.0);
    test_value(model3.inverse_cutoff(), 1.0e-3);
    test_value(model3.cutoff().GeV(), 1.0);

    // Test XML constructor
    GXml         xml(m_xml_model_point_einvplaw);
    GXmlElement* element = xml.element(0)->element(0)->element("spectrum", 0);
    GModelSpectralExpInvPlaw model4(*element);
    test_value(model4.size(), 4);
    test_value(model4.type(), "ExponentialCutoffPowerLaw", "Check type of model");
    test_value(model4.prefactor(), 5.7e-16);
    test_value(model4.index(), -2.48);
    test_value(model4.pivot().TeV(), 0.3);
    test_value(model4.inverse_cutoff(), 1.0e-6);
    test_value(model4.cutoff().TeV(), 1.0);

    // Test prefactor method
    model4.prefactor(2.3e-16);
    test_value(model4.prefactor(), 2.3e-16);

    // Test index method
    model4.index(-2.6);
    test_value(model4.index(), -2.6);

    // Test pivot method
    model4.pivot(GEnergy(0.5, "TeV"));
    test_value(model4.pivot().TeV(), 0.5);

    // Test cutoff parameter method
    model4.inverse_cutoff(4.2e-9);
    test_value(model4.inverse_cutoff(), 4.2e-9);

    // Test cutoff method
    model4.cutoff(GEnergy(10.0, "TeV"));
    test_value(model4.cutoff().TeV(), 10.0);

    // Test operator access
    const char* strarray[] = {"Prefactor", "Index", "PivotEnergy", "InverseCutoffEnergy"};
    for (int i = 0; i < 4; ++i) {
        std::string keyname(strarray[i]);
        model4[keyname].remove_range(); // To allow setting of any value
        model4[keyname].value(2.1);
        model4[keyname].error(1.9);
        model4[keyname].gradient(0.8);
        test_value(model4[keyname].value(), 2.1);
        test_value(model4[keyname].error(), 1.9);
        test_value(model4[keyname].gradient(), 0.8);
    }

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Test GModelSpectralSuperExpPlaw class
 ***************************************************************************/
void TestGModel::test_supeplaw(void)
{
    // Test void constructor
    GModelSpectralSuperExpPlaw model1;
    test_value(model1.type(), "SuperExponentialCutoffPowerLaw", "Check void mode type");

    // Test value constructor
    GModelSpectralSuperExpPlaw model2(2.0, -2.1, GEnergy(100.0, "MeV"),
                                      GEnergy(1.0, "GeV"), 1.1);
    test_value(model2.prefactor(), 2.0);
    test_value(model2.index1(), -2.1);
    test_value(model2.pivot().MeV(), 100.0);
    test_value(model2.cutoff().GeV(), 1.0);
    test_value(model2.index2(), 1.1);
    
    // Test XML constructor
    GXml                       xml(m_xml_model_point_supeplaw);
    GXmlElement*               element = xml.element(0)->element(0)->element("spectrum", 0);
    GModelSpectralSuperExpPlaw model3(*element);
    test_value(model3.size(), 5);
    test_value(model3.type(), "SuperExponentialCutoffPowerLaw", "Check model type");
    test_value(model3.prefactor(), 1e-16);
    test_value(model3.index1(), -2.0);
    test_value(model3.pivot().TeV(), 1.0);
    test_value(model3.cutoff().TeV(), 1.0);
    test_value(model3.index2(), 1.5);

    // Test prefactor method
    model3.prefactor(2.3e-16);
    test_value(model3.prefactor(), 2.3e-16);

    // Test index1 method
    model3.index1(-2.6);
    test_value(model3.index1(), -2.6);

    // Test pivot method
    model3.pivot(GEnergy(0.5, "TeV"));
    test_value(model3.pivot().TeV(), 0.5);

    // Test cutoff method
    model3.cutoff(GEnergy(10.0, "TeV"));
    test_value(model3.cutoff().TeV(), 10.0);

    // Test index2 method
    model3.index2(1.7);
    test_value(model3.index2(), 1.7);

    // Test operator access
    const char* strarray[] = {"Prefactor", "Index1", "PivotEnergy",
                              "CutoffEnergy", "Index2"};
    for (int i = 0; i < 5; ++i) {
        std::string keyname(strarray[i]);
        model3[keyname].remove_range(); // To allow setting of any value
        model3[keyname].value(2.1);
        model3[keyname].error(1.9);
        model3[keyname].gradient(0.8);
        test_value(model3[keyname].value(), 2.1);
        test_value(model3[keyname].error(), 1.9);
        test_value(model3[keyname].gradient(), 0.8);
    }

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Test GModelSpectralBrokenPlaw class
 ***************************************************************************/
void TestGModel::test_bplaw(void)
{
    // Test void constructor
    GModelSpectralBrokenPlaw model1;
    test_value(model1.type(), "BrokenPowerLaw", "Check void model type");

    // Test value constructor
    GModelSpectralBrokenPlaw model2(2.0, -2.1, GEnergy(100.0, "MeV"), -2.8);
    test_value(model2.prefactor(), 2.0);
    test_value(model2.index1(), -2.1);
    test_value(model2.breakenergy().MeV(), 100.0);
    test_value(model2.index2(), -2.8);
    
    // Test XML constructor
    GXml                     xml(m_xml_model_point_bplaw);
    GXmlElement*             element = xml.element(0)->element(0)->element("spectrum", 0);
    GModelSpectralBrokenPlaw model3(*element);
    test_value(model3.size(), 4);
    test_value(model3.type(), "BrokenPowerLaw", "Check model type");
    test_value(model3.prefactor(), 5.7e-16);
    test_value(model3.index1(), -2.48);
    test_value(model3.breakenergy().TeV(), 0.3);
    test_value(model3.index2(), -2.70);

    // Test prefactor method
    model3.prefactor(2.3e-16);
    test_value(model3.prefactor(), 2.3e-16);

    // Test index1 method
    model3.index1(-2.6);
    test_value(model3.index1(), -2.6);

    // Test breakenergy method
    model3.breakenergy(GEnergy(0.5, "TeV"));
    test_value(model3.breakenergy().TeV(), 0.5);

    // Test index2 method
    model3.index2(-3.6);
    test_value(model3.index2(), -3.6);

    // Test operator access
    const char* strarray[] = {"Prefactor", "Index1", "BreakEnergy", "Index2"};
    for (int i = 0; i < 4; ++i) {
        std::string keyname(strarray[i]);
        model3[keyname].remove_range(); // To allow setting of any value
        model3[keyname].value(2.1);
        model3[keyname].error(1.9);
        model3[keyname].gradient(0.8);
        test_value(model3[keyname].value(), 2.1);
        test_value(model3[keyname].error(), 1.9);
        test_value(model3[keyname].gradient(), 0.8);
    }

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Test GModelSpectralLogParabola class
 ***************************************************************************/
void TestGModel::test_logparabola(void)
{
    // Test void constructor
    GModelSpectralLogParabola model1;
    test_value(model1.type(), "LogParabola", "Check void model type");

    // Test value constructor
    GModelSpectralLogParabola model2(2.0, -2.1, GEnergy(100.0, "MeV"), -0.2);
    test_value(model2.prefactor(), 2.0);
    test_value(model2.index(), -2.1);
    test_value(model2.pivot().MeV(), 100.0);
    test_value(model2.curvature(), -0.2);
    
    // Test XML constructor
    GXml                      xml(m_xml_model_point_logparabola);
    GXmlElement*              element = xml.element(0)->element(0)->element("spectrum", 0);
    GModelSpectralLogParabola model3(*element);
    test_value(model3.size(), 4);
    test_value(model3.type(), "LogParabola", "Check model type");
    test_value(model3.prefactor(), 5.878e-17);
    test_value(model3.index(), -2.32473);
    test_value(model3.pivot().TeV(), 1.0);
    test_value(model3.curvature(), -0.074);

    // Test prefactor method
    model3.prefactor(2.3e-16);
    test_value(model3.prefactor(), 2.3e-16);

    // Test index method
    model3.index(-2.6);
    test_value(model3.index(), -2.6);

    // Test pivot method
    model3.pivot(GEnergy(0.5, "TeV"));
    test_value(model3.pivot().TeV(), 0.5);

    // Test curvature method
    model3.curvature(-0.1);
    test_value(model3.curvature(), -0.1);

    // Test operator access
    const char* strarray[] = {"Prefactor", "Index", "PivotEnergy", "Curvature"};
    for (int i = 0; i < 4; ++i) {
        std::string keyname(strarray[i]);
        model3[keyname].remove_range(); // To allow setting of any value
        model3[keyname].value(2.1);
        model3[keyname].error(1.9);
        model3[keyname].gradient(0.8);
        test_value(model3[keyname].value(), 2.1);
        test_value(model3[keyname].error(), 1.9);
        test_value(model3[keyname].gradient(), 0.8);
    }

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Test GModelSpectralNodes class
 ***************************************************************************/
void TestGModel::test_nodes(void)
{
    // Test void constructor
    GModelSpectralNodes model1;
    test_value(model1.type(), "NodeFunction", "Check void model type");

    // Test node function manimulation
    GModelSpectralNodes model2;
    model2.reserve(3);
    test_value(model2.size(), 0);
    test_value(model2.nodes(), 0);
    model2.append(GEnergy(1.0, "MeV"), 1.0);
    test_value(model2.size(), 2);
    test_value(model2.nodes(), 1);
    test_value(model2.energy(0).MeV(), 1.0);
    test_value(model2.intensity(0), 1.0);
    model2.append(GEnergy(10.0, "MeV"), 0.1);
    test_value(model2.size(), 4);
    test_value(model2.nodes(), 2);
    test_value(model2.energy(0).MeV(), 1.0);
    test_value(model2.energy(1).MeV(), 10.0);
    test_value(model2.intensity(0), 1.0);
    test_value(model2.intensity(1), 0.1);
    model2.remove(0);
    test_value(model2.size(), 2);
    test_value(model2.nodes(), 1);
    test_value(model2.energy(0).MeV(), 10.0);
    test_value(model2.intensity(0), 0.1);
    model2.insert(0, GEnergy(1.0, "MeV"), 1.0);
    test_value(model2.size(), 4);
    test_value(model2.nodes(), 2);
    test_value(model2.energy(0).MeV(), 1.0);
    test_value(model2.energy(1).MeV(), 10.0);
    test_value(model2.intensity(0), 1.0);
    test_value(model2.intensity(1), 0.1);
    model2.extend(model2);
    test_value(model2.size(), 8);
    test_value(model2.nodes(), 4);
    test_value(model2.energy(0).MeV(), 1.0);
    test_value(model2.energy(1).MeV(), 10.0);
    test_value(model2.energy(2).MeV(), 1.0);
    test_value(model2.energy(3).MeV(), 10.0);
    test_value(model2.intensity(0), 1.0);
    test_value(model2.intensity(1), 0.1);
    test_value(model2.intensity(2), 1.0);
    test_value(model2.intensity(3), 0.1);
   
    // Test XML constructor
    GXml                      xml(m_xml_model_point_nodes);
    GXmlElement*              element = xml.element(0)->element(0)->element("spectrum", 0);
    GModelSpectralNodes model3(*element);
    test_value(model3.size(), 4);
    test_value(model3.nodes(), 2);
    test_value(model3.type(), "NodeFunction", "Check model type");
    test_value(model3.energy(0).MeV(), 1.0);
    test_value(model3.energy(1).MeV(), 10.0);
    test_value(model3.intensity(0), 1.0e-7);
    test_value(model3.intensity(1), 0.1e-7);

    // Test energy method
    model3.energy(0, GEnergy(0.1, "MeV"));
    test_value(model3.energy(0).MeV(), 0.1);

    // Test intensity method
    model3.intensity(0, 2.0e-7);
    test_value(model3.intensity(0), 2.0e-7);

    // Test operator access
    const char* strarray[] = {"Energy0", "Energy1", "Intensity0", "Intensity1"};
    for (int i = 0; i < 4; ++i) {
        std::string keyname(strarray[i]);
        model3[keyname].remove_range(); // To allow setting of any value
        model3[keyname].value(2.1);
        model3[keyname].error(1.9);
        model3[keyname].gradient(0.8);
        test_value(model3[keyname].value(), 2.1);
        test_value(model3[keyname].error(), 1.9);
        test_value(model3[keyname].gradient(), 0.8);
    }

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Test GModelSpectralFunc class
 ***************************************************************************/
void TestGModel::test_filefct(void)
{
    // Test void constructor
    GModelSpectralFunc model1;
    test_value(model1.type(), "FileFunction", "Check void model type");

    // Test value constructor
    GModelSpectralFunc model2(m_filefct, 2.0);
    test_value(model2.filename().url(), m_filefct,
                "Check file function data file name");
    test_value(model2.norm(), 2.0);
   
    // Test XML constructor
    GXml               xml(m_xml_model_point_filefct);
    GXmlElement*       element = xml.element(0)->element(0)->element("spectrum", 0);
    GModelSpectralFunc model3(*element);
    test_value(model3.size(), 1);
    test_value(model3.type(), "FileFunction", "Check model type");
    test_value(model3.filename().url(), m_filefct,
               "Check file function data file name");
    test_value(model3.norm(), 1.0);

    // Test filename method
    model3.filename(m_filefct);
    test_value(model3.filename().url(), m_filefct,
               "Check file function data file name");

    // Test norm method
    model3.norm(3.0);
    test_value(model3.norm(), 3.0);

    // Test operator access
    const char* strarray[] = {"Normalization"};
    for (int i = 0; i < 1; ++i) {
        std::string keyname(strarray[i]);
        model3[keyname].remove_range(); // To allow setting of any value
        model3[keyname].value(2.1);
        model3[keyname].error(1.9);
        model3[keyname].gradient(0.8);
        test_value(model3[keyname].value(), 2.1);
        test_value(model3[keyname].error(), 1.9);
        test_value(model3[keyname].gradient(), 0.8);
    }

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Test GModelTemporalConst class
 ***************************************************************************/
void TestGModel::test_temp_const(void)
{
    // Test void constructor
    GModelTemporalConst model1;
    test_value(model1.type(), "Constant", "Check void model type");

    // Test value constructor
    GModelTemporalConst model2(3.0);
    test_value(model2.norm(), 3.0);
    
    // Test XML constructor
    GModelTemporalConst model3;

    // Test value method
    model3.norm(3.9);
    test_value(model3.norm(), 3.9);

    // Test operator access
    const char* strarray[] = {"Normalization"};
    for (int i = 0; i < 1; ++i) {
        std::string keyname(strarray[i]);
        model3[keyname].remove_range(); // To allow setting of any value
        model3[keyname].value(2.1);
        model3[keyname].error(1.9);
        model3[keyname].gradient(0.8);
        test_value(model3[keyname].value(), 2.1);
        test_value(model3[keyname].error(), 1.9);
        test_value(model3[keyname].gradient(), 0.8);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test loading and saving of XML model.
 *
 * @param[in] name Model name.
 * @param[in] filename XML model filename.
 *
 * Test loading and saving of XML model.
 ***************************************************************************/
void TestGModel::test_xml_model(const std::string& name,
                                const std::string& filename)
{
    // Set filename of test file
    std::string fname = "test_xml_" + name + ".xml";

    // Test saving and reloading
    GModels models(filename);
    models.save(fname);
    models.load(fname);
    models.save(fname);
    models.load(fname);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test spatial model XML reading and writing
 ***************************************************************************/
void TestGModel::test_spatial_model(void)
{
    // Test spatial models XML interface
    test_xml_model("GModelSpatialPointSource",     m_xml_model_point_plaw);
    test_xml_model("GModelSpatialRadialDisk",      m_xml_model_radial_disk);
    test_xml_model("GModelSpatialRadialGauss",     m_xml_model_radial_gauss);
    test_xml_model("GModelSpatialRadialShell",     m_xml_model_radial_shell);
    test_xml_model("GModelSpatialEllipticalDisk",  m_xml_model_elliptical_disk);
    test_xml_model("GModelSpatialEllipticalGauss", m_xml_model_elliptical_gauss);
    test_xml_model("GModelSpatialDiffuseConst",    m_xml_model_diffuse_const);
    test_xml_model("GModelSpatialDiffuseMap",      m_xml_model_diffuse_map);
    test_xml_model("GModelSpatialDiffuseCube",     m_xml_model_diffuse_cube);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test spectral model XML reading and writing
 ***************************************************************************/
void TestGModel::test_spectral_model(void)
{
    // Test spectral models XML interface
    test_xml_model("GModelSpectralConst",          m_xml_model_point_const);
    test_xml_model("GModelSpectralPlaw",           m_xml_model_point_plaw);
    test_xml_model("GModelSpectralPlawPhotonFlux", m_xml_model_point_plaw_phflux);
    test_xml_model("GModelSpectralPlawEnergyFlux", m_xml_model_point_plaw_eflux);
    test_xml_model("GModelSpectralExpPaw",         m_xml_model_point_eplaw);
    test_xml_model("GModelSpectralExpInvPaw",      m_xml_model_point_einvplaw);
    test_xml_model("GModelSpectralBrokenPlaw",     m_xml_model_point_bplaw);
    test_xml_model("GModelSpectralSuperExpPlaw",   m_xml_model_point_supeplaw);
    test_xml_model("GModelSpectralLogParabola",    m_xml_model_point_logparabola);
    test_xml_model("GModelSpectralNodes",          m_xml_model_point_nodes);
    test_xml_model("GModelSpectralFunc",           m_xml_model_point_filefct);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test model base class functions
 ***************************************************************************/
void TestGModel::test_model(void)
{
    // Setup some models (use the sky model for that)
    GModelSky model;

    // Test instruments setting
    model.instruments("CTA");
    test_assert(model.instruments() == "CTA",
                "Expected \"CTA\" instruments, found \""+
                model.instruments()+"\".");
    test_assert(model.is_valid("CTA", ""), "Model is not a \"CTA\" model.");
    test_assert(model.is_valid("CTA", "123"), "Model ID selection error.");
    model.instruments("CTA,LAT");
    test_assert(model.instruments() == "CTA,LAT",
                "Expected \"CTA,LAT\" instruments, found \""+
                model.instruments()+"\".");
    test_assert(model.is_valid("CTA", ""), "Model is not a \"CTA\" model.");
    test_assert(model.is_valid("LAT", ""), "Model is not a \"LAT\" model.");
    test_assert(model.is_valid("CTA", "123"), "Model ID selection error.");
    test_assert(model.is_valid("LAT", "123"), "Model ID selection error.");

    // Test ID setting
    model.clear();
    model.ids("123");
    test_assert(model.ids() == "123",
                "Expected ID \"123\" instruments, found \""+model.ids()+"\".");
    test_assert(model.is_valid("", "123"), "Model ID is not \"123\".");
    test_assert(model.is_valid("CTA", "123"), "Model instrument selection error.");
    model.ids("123,456");
    test_assert(model.ids() == "123,456",
                "Expected ID \"123,456\" instruments, found \""+
                model.ids()+"\".");
    test_assert(model.is_valid("", "123"), "Model ID is not \"123\".");
    test_assert(model.is_valid("", "456"), "Model ID is not \"456\".");
    test_assert(model.is_valid("CTA", "123"), "Model instrument selection error.");
    test_assert(model.is_valid("CTA", "456"), "Model instrument selection error.");

    // Test instrument and ID setting
    model.clear();
    model.instruments("CTA");
    model.ids("123");
    test_assert(model.instruments() == "CTA",
                "Expected \"CTA\" instruments, found \""+
                model.instruments()+"\".");
    test_assert(model.ids() == "123",
                "Expected ID \"123\" instruments, found \""+model.ids()+"\".");
    test_assert(model.is_valid("CTA", "123"), "Test model instrument and ID.");
    test_assert(!model.is_valid("CTAx", "123"), "Test failure mode.");
    test_assert(!model.is_valid("CTA", "123x"), "Test failure mode.");
    model.instruments("CTA,LAT");
    model.ids("123,456");
    test_assert(model.is_valid("", ""), "Test model instrument and ID.");
    test_assert(model.is_valid("CTA", "123"), "Test model instrument and ID.");
    test_assert(model.is_valid("LAT", "123"), "Test model instrument and ID.");
    test_assert(model.is_valid("CTA", "456"), "Test model instrument and ID.");
    test_assert(model.is_valid("LAT", "456"), "Test model instrument and ID.");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test model container handling
 ***************************************************************************/
void TestGModel::test_models(void)
{
    // Test void constructor
    test_try("Test void constructor");
    try {
        GModels models;
        test_assert(models.is_empty(), "Model container not empty.");
        test_value(models.size(), 0);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test load constructor
    test_try("Test load constructor");
    try {
        GModels models(m_xml_file);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test saving and loading
    test_try("Test saving and loading");
    try {
        GModels models(m_xml_file);
        models.save("test.xml");
        models.load("test.xml");
        models.save("test.xml");
        models.load("test.xml");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test model manipulation
    test_try("Test model access");
    try {
        GModels models(m_xml_file);
        GModel* model = models[0]->clone();
        models.reserve(5);
        test_assert(models.contains("1FGL J0005.7+3815"), 
                    "Model \"1FGL J0005.7+3815\" not found.");
        test_assert(!models.contains("2FGL J0005.7+3815"), 
                    "Model \"2FGL J0005.7+3815\" found but not expected.");
        test_assert(!models.is_empty(), "Model container is empty.");
        test_value(models.size(), 1);
        
        // Append model with same name
        test_try("Append model with same name");
        try {
            models.append(*(models[0]));
            test_try_failure("Appending of model with same name is not refused.");
        }
        catch (GException::invalid_value &e) {
            test_try_success();
        }
        catch (std::exception &e) {
            test_try_failure(e);
        }
        
        // Append model with different name
        model->name("Appended model");
        models.append(*model);

        // Insert model with same name
        test_try("Insert model with same name");
        try {
            models.insert(0, *(models.at(0)));
            test_try_failure("Inserting of model with same name is not refused.");
        }
        catch (GException::invalid_value &e) {
            test_try_success();
        }
        catch (std::exception &e) {
            test_try_failure(e);
        }

        // Insert model with different name
        model->name("Inserted model");
        models.insert(0, *model);

        // Set model with same name
        test_try("Set model with same name");
        try {
            models.set(0, *(models.at(1)));
            test_try_failure("Setting of model with same name is not refused.");
        }
        catch (GException::invalid_value &e) {
            test_try_success();
        }
        catch (std::exception &e) {
            test_try_failure(e);
        }

        // Set model with different name
        model->name("Set model");
        models.set(0, *model);

        // Make sure that we have now 3 models in container
        test_assert(!models.is_empty(), "Model container is empty.");
        test_value(models.size(), 3);

        // Remove "Set model"
        models.remove(0);
        test_assert(!models.is_empty(), "Model container is empty.");
        test_value(models.size(), 2);

        // Remove "1FGL J0005.7+3815"
        models.remove(0);
        test_assert(!models.is_empty(), "Model container is empty.");
        test_value(models.size(), 1);

        // Remove "Appended model"
        models.remove("Appended model");
        test_assert(models.is_empty(), "Model container is not empty.");
        test_value(models.size(), 0);

        // Reload model
        models.load(m_xml_file);
        test_assert(!models.is_empty(), "Model container is empty.");
        test_value(models.size(), 1);
        
        // Append identical container
        test_try("Append identical container");
        try {
            models.extend(models);
            test_try_failure("Appending of identical container is not refused.");
        }
        catch (GException::invalid_value &e) {
            test_try_success();
        }
        catch (std::exception &e) {
            test_try_failure(e);
        }

        // Append different container
        GModels other_models = models;
        other_models[0]->name("New name");
        models.extend(other_models);
        test_assert(!models.is_empty(), "Model container is empty.");
        test_value(models.size(), 2);

        // Free model
        delete model;

        // Signal success
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }
    

    // Setup Crab model
    GModelSky crab;
    test_try("Setup Crab model");
    try {
        GModelSpatialPointSource point_source(83.6331, +22.0145);
        GModelSpectralPlaw power_law(2.0, -2.1, GEnergy(100.0, "MeV"));
        crab = GModelSky(point_source, power_law);
        crab.name("Crab");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Set model scaling
    GModelPar lat("LAT", 1.0);
    GModelPar cta("CTA", 0.5);
    crab.scale(lat);
    crab.scale(lat); // In purpose to check if parameter is appended only once
    crab.scale(cta);

    // Test model scaling
    test_value(crab.scale("LAT").value(), 1.0);
    test_value(crab.scale("CTA").value(), 0.5);
    test_value(crab.scale("COM").value(), 1.0);

    // Test saving and loading
    GModels models;
    models.append(crab);
    models.save("test_instrument.xml");
    models.load("test_instrument.xml");
    test_value(models[0]->scale("LAT").value(), 1.0);
    test_value(models[0]->scale("CTA").value(), 0.5);
    test_value(models[0]->scale("COM").value(), 1.0);

    // Exit test
    return;

}


/***********************************************************************//**
 * @brief Test model registry
 ***************************************************************************/
void TestGModel::test_model_registry(void)
{
    // Test GModelRegistry allocators
    test_try("Test GModelRegistry model allocation");
    try {
        // Loop over all models in registry
        GModelRegistry registry;
        int num = registry.size();
        for (int i = 0; i < num; ++i) {
            GModel* ptr = registry.alloc(registry.name(i));
            test_assert(ptr != NULL, "Model pointer for \"" +
                                     registry.name(i)+"\" is NULL");
            if (ptr != NULL) {
            
                // Test model type
                test_assert(ptr->type() == registry.name(i),
                            "Expected \""+registry.name(i)+"\" instead"
                            " of \""+ptr->type()+"\".");

                // Free model
                delete ptr;
            }
        } // endfor: looped over all models
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test GModelSpatialRegistry allocator for invalid XML element
    test_try("Test GModelSpatialRegistry allocator for invalid XML element");
    try {
        GXmlElement           xml;
        GModelSpatialRegistry registry;
        GModelSpatial*        ptr = registry.alloc(xml);
        test_try_failure("Invalid XML element shall throw an exception");
    }
    catch (GException::invalid_value &e) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test GModelSpectralRegistry allocator for invalid XML element
    test_try("Test GModelSpectralRegistry allocator for invalid XML element");
    try {
        GXmlElement            xml;
        GModelSpectralRegistry registry;
        GModelSpectral*        ptr = registry.alloc(xml);
        test_try_failure("Invalid XML element shall throw an exception");
    }
    catch (GException::invalid_value &e) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test GModelTemporalRegistry allocator for invalid XML element
    test_try("Test GModelTemporalRegistry allocator for invalid XML element");
    try {
        GXmlElement            xml;
        GModelTemporalRegistry registry;
        GModelTemporal*        ptr = registry.alloc(xml);
        test_try_failure("Invalid XML element shall throw an exception");
    }
    catch (GException::invalid_value &e) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Exit test
    return;

}


/***********************************************************************//**
 * @brief Test radial disk legacy model
 ***************************************************************************/
void TestGModel::test_legacy_model_radial_disk(void)
{
    // Load model from XML file
    GModels models(m_xml_legacy_radial_disk);

    // Extract spatial component
    GModelSpatialRadialDisk* spatial = static_cast<GModelSpatialRadialDisk*>(
                                         static_cast<GModelSky*>(
                                           models[0])->spatial());

    // Test model values
    test_value(spatial->type(), "DiskFunction", "Check model type");
    test_value(spatial->ra(), 83.6331, 1.0e-7, "Check Right Ascension");
    test_value(spatial->dec(), 22.0145, 1.0e-7, "Check Declination");
    test_value(spatial->radius(), 0.45, 1.0e-7, "Check radius");

    // Save file to disk and reload file (tests the proper saving and loading)
    models.save("test_xml_legacy_radial_disk.xml");
    models.clear();
    models.load("test_xml_legacy_radial_disk.xml");

    // Extract spatial component
    spatial = static_cast<GModelSpatialRadialDisk*>(
                static_cast<GModelSky*>(models[0])->spatial());

    // Re-test model values
    test_value(spatial->type(), "DiskFunction", "Check model type");
    test_value(spatial->ra(), 83.6331, 1.0e-7, "Check Right Ascension");
    test_value(spatial->dec(), 22.0145, 1.0e-7, "Check Declination");
    test_value(spatial->radius(), 0.45, 1.0e-7, "Check radius");

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Test radial Gaussian legacy model
 ***************************************************************************/
void TestGModel::test_legacy_model_radial_gauss(void)
{
    // Load model from XML file
    GModels models(m_xml_legacy_radial_gauss);

    // Extract spatial component
    GModelSpatialRadialGauss* spatial = static_cast<GModelSpatialRadialGauss*>(
                                          static_cast<GModelSky*>(
                                            models[0])->spatial());

    // Test model values
    test_value(spatial->type(), "GaussFunction", "Check model type");
    test_value(spatial->ra(), 83.6331, 1.0e-7, "Check Right Ascension");
    test_value(spatial->dec(), 22.0145, 1.0e-7, "Check Declination");
    test_value(spatial->sigma(), 0.20, 1.0e-7, "Check Sigma");

    // Save file to disk and reload file (tests the proper saving and loading)
    models.save("test_xml_legacy_radial_gauss.xml");
    models.clear();
    models.load("test_xml_legacy_radial_gauss.xml");

    // Extract spatial component
    spatial = static_cast<GModelSpatialRadialGauss*>(
                static_cast<GModelSky*>(models[0])->spatial());

    // Re-test model values
    test_value(spatial->type(), "GaussFunction", "Check model type");
    test_value(spatial->ra(), 83.6331, 1.0e-7, "Check Right Ascension");
    test_value(spatial->dec(), 22.0145, 1.0e-7, "Check Declination");
    test_value(spatial->sigma(), 0.20, 1.0e-7, "Check Sigma");

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Test radial shell legacy model
 ***************************************************************************/
void TestGModel::test_legacy_model_radial_shell(void)
{
    // Load model from XML file
    GModels models(m_xml_legacy_radial_shell);

    // Extract spatial component
    GModelSpatialRadialShell* spatial = static_cast<GModelSpatialRadialShell*>(
                                          static_cast<GModelSky*>(
                                            models[0])->spatial());

    // Test model values
    test_value(spatial->type(), "ShellFunction", "Check model type");
    test_value(spatial->ra(), 83.6331, 1.0e-7, "Check Right Ascension");
    test_value(spatial->dec(), 22.0145, 1.0e-7, "Check Declination");
    test_value(spatial->radius(), 0.30, 1.0e-7, "Check radius");
    test_value(spatial->width(), 0.10, 1.0e-7, "Check width");

    // Save file to disk and reload file (tests the proper saving and loading)
    models.save("test_xml_legacy_radial_shell.xml");
    models.clear();
    models.load("test_xml_legacy_radial_shell.xml");

    // Extract spatial component
    spatial = static_cast<GModelSpatialRadialShell*>(
                static_cast<GModelSky*>(models[0])->spatial());

    // Re-test model values
    test_value(spatial->type(), "ShellFunction", "Check model type");
    test_value(spatial->ra(), 83.6331, 1.0e-7, "Check Right Ascension");
    test_value(spatial->dec(), 22.0145, 1.0e-7, "Check Declination");
    test_value(spatial->radius(), 0.30, 1.0e-7, "Check radius");
    test_value(spatial->width(), 0.10, 1.0e-7, "Check width");

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Test elliptical Gaussian legacy model
 ***************************************************************************/
void TestGModel::test_legacy_model_elliptical_gauss(void)
{
    // Load model from XML file
    GModels models(m_xml_legacy_elliptical_gauss);

    // Extract spatial component
    GModelSpatialEllipticalGauss* spatial = static_cast<GModelSpatialEllipticalGauss*>(
                                              static_cast<GModelSky*>(
                                                models[0])->spatial());

    // Test model values
    test_value(spatial->type(), "EllipticalGauss", "Check model type");
    test_value(spatial->ra(), 83.6331, 1.0e-7, "Check Right Ascension");
    test_value(spatial->dec(), 22.0145, 1.0e-7, "Check Declination");
    test_value(spatial->posangle(), 45.0, 1.0e-7, "Check position angle");
    test_value(spatial->semimajor(), 0.3, 1.0e-7, "Check semi major axis");
    test_value(spatial->semiminor(), 0.1, 1.0e-7, "Check semi minor axis");

    // Save file to disk and reload file (tests the proper saving and loading)
    models.save("test_xml_legacy_elliptical_gauss.xml");
    models.clear();
    models.load("test_xml_legacy_elliptical_gauss.xml");

    // Extract spatial component
    spatial = static_cast<GModelSpatialEllipticalGauss*>(
                static_cast<GModelSky*>(models[0])->spatial());

    // Re-test model values
    test_value(spatial->type(), "EllipticalGauss", "Check model type");
    test_value(spatial->ra(), 83.6331, 1.0e-7, "Check Right Ascension");
    test_value(spatial->dec(), 22.0145, 1.0e-7, "Check Declination");
    test_value(spatial->posangle(), 45.0, 1.0e-7, "Check position angle");
    test_value(spatial->semimajor(), 0.3, 1.0e-7, "Check semi major axis");
    test_value(spatial->semiminor(), 0.1, 1.0e-7, "Check semi minor axis");

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Test diffuse constant legacy model
 ***************************************************************************/
void TestGModel::test_legacy_model_diffuse_const(void)
{
    // Load model from XML file
    GModels models(m_xml_legacy_diffuse_const);

    // Extract spatial component
    GModelSpatialDiffuseConst* spatial = static_cast<GModelSpatialDiffuseConst*>(
                                           static_cast<GModelSky*>(
                                             models[0])->spatial());

    // Test model values
    test_value(spatial->type(), "ConstantValue", "Check model type");
    test_value(spatial->value(), 1.0, 1.0e-7, "Check constant value");

    // Save file to disk and reload file (tests the proper saving and loading)
    models.save("test_xml_legacy_diffuse_const.xml");
    models.clear();
    models.load("test_xml_legacy_diffuse_const.xml");

    // Extract spatial component
    spatial = static_cast<GModelSpatialDiffuseConst*>(
                static_cast<GModelSky*>(models[0])->spatial());

    // Re-test model values
    test_value(spatial->type(), "ConstantValue", "Check model type");
    test_value(spatial->value(), 1.0, 1.0e-7, "Check constant value");

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Test diffuse map legacy model
 ***************************************************************************/
void TestGModel::test_legacy_model_diffuse_map(void)
{
    // Load model from XML file
    GModels models(m_xml_legacy_diffuse_map);

    // Extract spatial component
    GModelSpatialDiffuseMap* spatial = static_cast<GModelSpatialDiffuseMap*>(
                                         static_cast<GModelSky*>(
                                           models[0])->spatial());

    // Test model values
    test_value(spatial->type(), "SpatialMap", "Check model type");
    test_value(spatial->value(), 1.0, 1.0e-7, "Check constant value");

    // Save file to disk and reload file (tests the proper saving and loading)
    models.save("test_xml_legacy_diffuse_map.xml");
    models.clear();
    models.load("test_xml_legacy_diffuse_map.xml");

    // Extract spatial component
    spatial = static_cast<GModelSpatialDiffuseMap*>(
                static_cast<GModelSky*>(models[0])->spatial());

    // Re-test model values
    test_value(spatial->type(), "SpatialMap", "Check model type");
    test_value(spatial->value(), 1.0, 1.0e-7, "Check constant value");

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Test diffuse map cube legacy model
 ***************************************************************************/
void TestGModel::test_legacy_model_diffuse_cube(void)
{
    // Load model from XML file
    GModels models(m_xml_legacy_diffuse_cube);

    // Extract spatial component
    GModelSpatialDiffuseCube* spatial = static_cast<GModelSpatialDiffuseCube*>(
                                          static_cast<GModelSky*>(
                                            models[0])->spatial());

    // Test model values
    test_value(spatial->type(), "MapCubeFunction", "Check model type");
    test_value(spatial->value(), 1.0, 1.0e-7, "Check constant value");

    // Save file to disk and reload file (tests the proper saving and loading)
    models.save("test_xml_legacy_diffuse_cube.xml");
    models.clear();
    models.load("test_xml_legacy_diffuse_cube.xml");

    // Extract spatial component
    spatial = static_cast<GModelSpatialDiffuseCube*>(
                static_cast<GModelSky*>(models[0])->spatial());

    // Re-test model values
    test_value(spatial->type(), "MapCubeFunction", "Check model type");
    test_value(spatial->value(), 1.0, 1.0e-7, "Check constant value");

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Test constant legacy model
 ***************************************************************************/
void TestGModel::test_legacy_model_point_const(void)
{
    // Load model from XML file
    GModels models(m_xml_legacy_point_const);

    // Extract spatial component
    GModelSpatialPointSource* spatial = static_cast<GModelSpatialPointSource*>(
                                          static_cast<GModelSky*>(
                                            models[0])->spatial());

    // Extract spectral component
    GModelSpectralConst* spectral = static_cast<GModelSpectralConst*>(
                                      static_cast<GModelSky*>(
                                        models[0])->spectral());

    // Test model values
    test_value(spatial->type(), "SkyDirFunction", "Check spatial model type");
    test_value(spatial->ra(), 83.6331, 1.0e-7, "Check Right Ascension");
    test_value(spatial->dec(), 22.0145, 1.0e-7, "Check Declination");
    test_value(spectral->type(), "ConstantValue", "Check spectral model type");
    test_value(spectral->value(), 5.7e-16, 1.0e-7, "Check constant value");

    // Save file to disk and reload file (tests the proper saving and loading)
    models.save("test_xml_legacy_point_const.xml");
    models.clear();
    models.load("test_xml_legacy_point_const.xml");

    // Extract spatial component
    spatial = static_cast<GModelSpatialPointSource*>(
                static_cast<GModelSky*>(models[0])->spatial());

    // Extract spectral component
    spectral = static_cast<GModelSpectralConst*>(
                 static_cast<GModelSky*>(models[0])->spectral());

    // Re-test model values
    test_value(spatial->type(), "SkyDirFunction", "Check spatial model type");
    test_value(spatial->ra(), 83.6331, 1.0e-7, "Check Right Ascension");
    test_value(spatial->dec(), 22.0145, 1.0e-7, "Check Declination");
    test_value(spectral->type(), "ConstantValue", "Check spectral model type");
    test_value(spectral->value(), 5.7e-16, 1.0e-7, "Check constant value");

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Test power law legacy model
 ***************************************************************************/
void TestGModel::test_legacy_model_point_plaw(void)
{
    // Load model from XML file
    GModels models(m_xml_legacy_point_plaw);

    // Extract spectral component
    GModelSpectralPlaw* spectral = static_cast<GModelSpectralPlaw*>(
                                     static_cast<GModelSky*>(
                                       models[0])->spectral());

    // Test model values
    test_value(spectral->type(), "PowerLaw", "Check model type");
    test_value(spectral->prefactor(), 5.7e-16, 1.0e-7, "Check prefactor");
    test_value(spectral->index(), -2.48, 1.0e-7, "Check index");
    test_value(spectral->pivot().GeV(), 300.0, 1.0e-7, "Check pivot energy");

    // Save file to disk and reload file (tests the proper saving and loading)
    models.save("test_xml_legacy_point_plaw.xml");
    models.clear();
    models.load("test_xml_legacy_point_plaw.xml");

    // Extract spectral component
    spectral = static_cast<GModelSpectralPlaw*>(
                 static_cast<GModelSky*>(models[0])->spectral());

    // Re-test model values
    test_value(spectral->type(), "PowerLaw", "Check model type");
    test_value(spectral->prefactor(), 5.7e-16, 1.0e-7, "Check prefactor");
    test_value(spectral->index(), -2.48, 1.0e-7, "Check index");
    test_value(spectral->pivot().GeV(), 300.0, 1.0e-7, "Check pivot energy");

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Test power law2 legacy model
 ***************************************************************************/
void TestGModel::test_legacy_model_point_plaw2(void)
{
    // Load model from XML file
    GModels models(m_xml_legacy_point_plaw2);

    // Extract spectral component
    GModelSpectralPlawPhotonFlux* spectral = static_cast<GModelSpectralPlawPhotonFlux*>(
                                      static_cast<GModelSky*>(
                                        models[0])->spectral());

    // Test model values
    test_value(spectral->type(), "PowerLaw2", "Check model type");
    test_value(spectral->flux(), 1.0e-7, 1.0e-7, "Check integral photon flux");
    test_value(spectral->index(), -2.0, 1.0e-7, "Check index");
    test_value(spectral->emin().MeV(), 100.0, 1.0e-7, "Check minimum energy");
    test_value(spectral->emax().MeV(), 500000.0, 1.0e-7, "Check maximum energy");

    // Save file to disk and reload file (tests the proper saving and loading)
    models.save("test_xml_legacy_point_plaw2.xml");
    models.clear();
    models.load("test_xml_legacy_point_plaw2.xml");

    // Extract spectral component
    spectral = static_cast<GModelSpectralPlawPhotonFlux*>(
                 static_cast<GModelSky*>(models[0])->spectral());

    // Re-test model values
    test_value(spectral->type(), "PowerLaw2", "Check model type");
    test_value(spectral->flux(), 1.0e-7, 1.0e-7, "Check integral photon flux");
    test_value(spectral->index(), -2.0, 1.0e-7, "Check index");
    test_value(spectral->emin().MeV(), 100.0, 1.0e-7, "Check minimum energy");
    test_value(spectral->emax().MeV(), 500000.0, 1.0e-7, "Check maximum energy");

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Test exponentially cut off power law legacy model
 ***************************************************************************/
void TestGModel::test_legacy_model_point_eplaw(void)
{
    // Load model from XML file
    GModels models(m_xml_legacy_point_eplaw);

    // Extract spectral component
    GModelSpectralExpPlaw* spectral = static_cast<GModelSpectralExpPlaw*>(
                                        static_cast<GModelSky*>(
                                          models[0])->spectral());

    // Test model values
    test_value(spectral->type(), "ExpCutoff", "Check model type");
    test_value(spectral->prefactor(), 5.7e-16, 1.0e-7, "Check prefactor");
    test_value(spectral->index(), -2.48, 1.0e-7, "Check index");
    test_value(spectral->pivot().GeV(), 300.0, 1.0e-7, "Check pivot energy");
    test_value(spectral->cutoff().TeV(), 1.0, 1.0e-7, "Check cut off energy");

    // Save file to disk and reload file (tests the proper saving and loading)
    models.save("test_xml_legacy_point_eplaw.xml");
    models.clear();
    models.load("test_xml_legacy_point_eplaw.xml");

    // Extract spectral component
    spectral = static_cast<GModelSpectralExpPlaw*>(
                 static_cast<GModelSky*>(models[0])->spectral());

    // Re-test model values
    test_value(spectral->type(), "ExpCutoff", "Check model type");
    test_value(spectral->prefactor(), 5.7e-16, 1.0e-7, "Check prefactor");
    test_value(spectral->index(), -2.48, 1.0e-7, "Check index");
    test_value(spectral->pivot().GeV(), 300.0, 1.0e-7, "Check pivot energy");
    test_value(spectral->cutoff().TeV(), 1.0, 1.0e-7, "Check cut off energy");

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Test super exponentially cut off power law legacy model
 ***************************************************************************/
void TestGModel::test_legacy_model_point_supeplaw(void)
{
    // Load model from XML file
    GModels models(m_xml_legacy_point_supeplaw);

    // Extract spectral component
    GModelSpectralSuperExpPlaw* spectral = static_cast<GModelSpectralSuperExpPlaw*>(
                                             static_cast<GModelSky*>(
                                               models[0])->spectral());

    // Test model values
    test_value(spectral->type(), "PLSuperExpCutoff", "Check model type");
    test_value(spectral->prefactor(), 1e-16, 1.0e-7, "Check prefactor");
    test_value(spectral->index1(), -2.0, 1.0e-7, "Check index 1");
    test_value(spectral->pivot().TeV(), 1.0, 1.0e-7, "Check pivot energy");
    test_value(spectral->index2(), 1.5, 1.0e-7, "Check index 2");
    test_value(spectral->cutoff().TeV(), 1.0, 1.0e-7, "Check cut off energy");

    // Save file to disk and reload file (tests the proper saving and loading)
    models.save("test_xml_legacy_point_supeplaw.xml");
    models.clear();
    models.load("test_xml_legacy_point_supeplaw.xml");

    // Extract spectral component
    spectral = static_cast<GModelSpectralSuperExpPlaw*>(
                 static_cast<GModelSky*>(models[0])->spectral());

    // Re-test model values
    test_value(spectral->type(), "PLSuperExpCutoff", "Check model type");
    test_value(spectral->prefactor(), 1e-16, 1.0e-7, "Check prefactor");
    test_value(spectral->index1(), -2.0, 1.0e-7, "Check index 1");
    test_value(spectral->pivot().TeV(), 1.0, 1.0e-7, "Check pivot energy");
    test_value(spectral->index2(), 1.5, 1.0e-7, "Check index 2");
    test_value(spectral->cutoff().TeV(), 1.0, 1.0e-7, "Check cut off energy");

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Test log parabola legacy model
 ***************************************************************************/
void TestGModel::test_legacy_model_point_logparabola(void)
{
    // Load model from XML file
    GModels models(m_xml_legacy_point_logparabola);

    // Extract spectral component
    GModelSpectralLogParabola* spectral = static_cast<GModelSpectralLogParabola*>(
                                            static_cast<GModelSky*>(
                                              models[0])->spectral());

    // Test model values
    test_value(spectral->type(), "LogParabola", "Check model type");
    test_value(spectral->prefactor(), 5.878e-16, 1.0e-7, "Check prefactor");
    test_value(spectral->index(), -2.32473, 1.0e-7, "Check index");
    test_value(spectral->pivot().TeV(), 1.0, 1.0e-7, "Check pivot energy");
    test_value(spectral->curvature(), -0.074, 1.0e-7, "Check curvature)");

    // Save file to disk and reload file (tests the proper saving and loading)
    models.save("test_xml_legacy_point_logparabola.xml");
    models.clear();
    models.load("test_xml_legacy_point_logparabola.xml");

    // Extract spectral component
    spectral = static_cast<GModelSpectralLogParabola*>(
                 static_cast<GModelSky*>(models[0])->spectral());

    // Re-test model values
    test_value(spectral->type(), "LogParabola", "Check model type");
    test_value(spectral->prefactor(), 5.878e-16, 1.0e-7, "Check prefactor");
    test_value(spectral->index(), -2.32473, 1.0e-7, "Check index");
    test_value(spectral->pivot().TeV(), 1.0, 1.0e-7, "Check pivot energy");
    test_value(spectral->curvature(), -0.074, 1.0e-7, "Check curvature)");

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Main test function.
 ***************************************************************************/
int main(void)
{
    // Create test suites
    GTestSuites testsuite("GModel");

    // Initialise success flag
    bool success = true;

    // Create a test suite
    TestGModel test;

    // Append test to the container
    testsuite.append(test);

    // Run the testsuites
    success = testsuite.run();

    // Save xml report
    testsuite.save("reports/GModel.xml");

    // Return
    return success ? 0:1;
}
