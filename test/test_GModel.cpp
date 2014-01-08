/***************************************************************************
 *                   test_GModel.cpp - test GModel class                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2013 by Juergen Knoedlseder                         *
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
#include <stdlib.h>
#include "test_GModel.hpp"
#include "GTools.hpp"


/***********************************************************************//**
 * @brief Set tests
 ***************************************************************************/
void TestGModel::set(void)
{
    // Set name
    name("GModel");

    // Set attributes
    m_map_file                    = "data/cena_lobes_parkes.fits";
    m_cube_file                   = "data/test_cube.fits";
    m_xml_file                    = "data/crab.xml";
    m_xml_model_point_const       = "data/model_point_const.xml";
    m_xml_model_point_plaw        = "data/model_point_plaw.xml";
    m_xml_model_point_plaw2       = "data/model_point_plaw2.xml";
    m_xml_model_point_eplaw       = "data/model_point_eplaw.xml";
    m_xml_model_point_bplaw       = "data/model_point_bplaw.xml";
    m_xml_model_point_logparabola = "data/model_point_logparabola.xml";
    m_xml_model_point_nodes       = "data/model_point_nodes.xml";
    m_xml_model_point_filefct     = "data/model_point_filefct.xml";
    m_xml_model_diffuse_const     = "data/model_diffuse_const.xml";
    m_xml_model_diffuse_cube      = "data/model_diffuse_cube.xml";
    m_xml_model_diffuse_map       = "data/model_diffuse_map.xml";
    m_xml_model_radial_disk       = "data/model_radial_disk.xml";
    m_xml_model_radial_gauss      = "data/model_radial_gauss.xml";
    m_xml_model_radial_shell      = "data/model_radial_shell.xml";
    m_xml_model_elliptical_disk   = "data/model_elliptical_disk.xml";

    // Append tests
    append(static_cast<pfunction>(&TestGModel::test_model_par), "Test GModelPar");
    append(static_cast<pfunction>(&TestGModel::test_sky_model), "Test GModelSky");

    // Append spatial model tests
    append(static_cast<pfunction>(&TestGModel::test_point_source), "Test GModelSpatialPointSource");
    append(static_cast<pfunction>(&TestGModel::test_radial_disk), "Test GModelSpatialRadialDisk");
    append(static_cast<pfunction>(&TestGModel::test_radial_gauss), "Test GModelSpatialRadialGauss");
    append(static_cast<pfunction>(&TestGModel::test_radial_shell), "Test GModelSpatialRadialShell");
    append(static_cast<pfunction>(&TestGModel::test_elliptical_disk), "Test GModelSpatialEllipticalDisk");
    append(static_cast<pfunction>(&TestGModel::test_diffuse_const), "Test GModelSpatialDiffuseConst");
    append(static_cast<pfunction>(&TestGModel::test_diffuse_cube), "Test GModelSpatialDiffuseCube");
    append(static_cast<pfunction>(&TestGModel::test_diffuse_map), "Test GModelSpatialDiffuseMap");
    append(static_cast<pfunction>(&TestGModel::test_spatial_model), "Test spatial model XML I/O");

    // Append spectral model tests
    append(static_cast<pfunction>(&TestGModel::test_const), "Test GModelSpectralConst");
    append(static_cast<pfunction>(&TestGModel::test_plaw), "Test GModelSpectralPlaw");
    append(static_cast<pfunction>(&TestGModel::test_plaw2), "Test GModelSpectralPlaw2");
    append(static_cast<pfunction>(&TestGModel::test_eplaw), "Test GModelSpectralExpPlaw");
    append(static_cast<pfunction>(&TestGModel::test_bplaw), "Test GModelSpectralBrokenPlaw");
    append(static_cast<pfunction>(&TestGModel::test_logparabola), "Test GModelSpectralLogParabola");
    append(static_cast<pfunction>(&TestGModel::test_nodes), "Test GModelSpectralNodes");
    append(static_cast<pfunction>(&TestGModel::test_filefct), "Test GModelSpectralFunc");
    append(static_cast<pfunction>(&TestGModel::test_spectral_model), "Test spectral model XML I/O");

    // Append temporal model tests
    append(static_cast<pfunction>(&TestGModel::test_temp_const), "Test GModelTemporalConst");

    // Append model container tests
    append(static_cast<pfunction>(&TestGModel::test_models), "Test GModels");

    // Append model registry tests
    append(static_cast<pfunction>(&TestGModel::test_model_registry), "Test model registries");

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
        test_value(par.gradient(), 102.0);
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

    // Test rescaling
    test_try("Test rescaling");
    try {
        GModelPar par("Test parameter", 3.0);
        par.error(3.0);
        par.gradient(3.0);
        par.min(3.0);
        par.max(3.0);
        par.scale(1.0);
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
        test_value(par.factor_gradient(), 1.0);
        test_value(par.factor_min(), 1.0);
        test_value(par.factor_max(), 1.0);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Exit test
    return;
}


/***********************************************************************//**
 * @brief Test GModelSky class
 ***************************************************************************/
void TestGModel::test_sky_model(void)
{
    // Test void constructor
    test_try("Test void constructor");
    try {
        GModelSky sky;
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test type constructor
    test_try("Test type constructor");
    try {
        GModelSky sky("My type");
        test_assert(sky.type() == "My type", "Model type \"My type\" expected.");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }
    
    // Test XML constructor, value and gradients
    test_try("Test XML constructor, value and gradients");
    try {
        // Test XML constructor
        GXml         xml(m_xml_file);
        GXmlElement* element = xml.element(0)->element(0);
        GModelSky    sky(*element);
        test_value(sky.size(), 6);
        test_assert(sky.name() == "1FGL J0005.7+3815",
                    "Expected source name \"1FGL J0005.7+3815\"");
        test_assert(sky.instruments() == "", "Expected no instruments");
        test_assert(sky.ids() == "", "Expected no observation identifiers");
        test_assert(sky.type() == "PointSource", "Expected \"PointSource\"");
        test_assert(sky.spatial() != NULL, "Expected spatial component");
        test_assert(sky.spectral() != NULL, "Expected spectral component");
        test_assert(sky.temporal() != NULL, "Expected temporal component");

        // Test value
        GSkyDir dir;
        dir.radec_deg(83.6331, 22.0145);
        GEnergy energy(100.0, "MeV");
        GTime   time(0.0);
        test_value(sky.value(GPhoton(dir, energy, time)), 1.73e-07);

        // Test gradient
        GVector vector = sky.gradients(GPhoton(dir, energy, time));
        test_value(vector[0], 0.0);
        test_value(vector[1], 0.0);
        test_value(vector[2], 1e-07);
        test_value(vector[3], 0.0);
        test_value(vector[4], 0.0);
        test_value(vector[5], 0.0);

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
 * @brief Test GModelSpatialPointSource class
 ***************************************************************************/
void TestGModel::test_point_source(void)
{
    // Test void constructor
    test_try("Test void constructor");
    try {
        GModelSpatialPointSource model;
        test_assert(model.type() == "SkyDirFunction",
                                    "Model type \"SkyDirFunction\" expected.");
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
        test_assert(model.type() == "SkyDirFunction", "Expected \"SkyDirFunction\"");
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
        test_assert(model.type() == "ConstantValue",
                                    "Model type \"ConstantValue\" expected.");
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
        test_assert(model.type() == "ConstantValue", "Expected \"ConstantValue\"");
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
        test_assert(model.type() == "MapCubeFunction",
                                    "Model type \"MapCubeFunction\" expected.");
        test_assert(model.filename() == "", "Model filename \"\" expected.");
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
        test_assert(model.filename() == m_cube_file, "Expected \""+m_cube_file+"\"");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test skymap value constructor
    test_try("Test skymap value constructor");
    try {
        GSkymap map("GAL", 16, "RING", 10);
        GEnergies energies;
        for (int i = 0; i < 10; ++i) {
            energies.append(GEnergy(double(i+1.0), "MeV"));
        }
        GModelSpatialDiffuseCube model(map, energies, 3.0);
        test_value(model.value(), 3.0);
        test_assert(model.filename() == "", "Expected \"\"");
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
        test_assert(model.type() == "MapCubeFunction", "Expected \"MapCubeFunction\"");
        test_value(model.value(), 1.0);
        test_assert(model.filename() == m_cube_file,
                                        "Model filename \""+m_cube_file+"\" expected.");

        // Test value method
        model.value(3.9);
        test_value(model.value(), 3.9);

        // Test filename method
        model.filename("Help me!");
        test_assert(model.filename() == "Help me!",
                                        "Model filename \"Help me!\" expected.");

        // Test cube method
        model.cube(GSkymap("GAL", 16, "RING", 10));
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
    test_try("Test void constructor");
    try {
        GModelSpatialDiffuseMap model;
        test_assert(model.type() == "SpatialMap",
                                    "Model type \"SpatialMap\" expected.");
        test_assert(model.filename() == "", "Model filename \"\" expected.");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test filename value constructor
    test_try("Test filename value constructor");
    try {
        GModelSpatialDiffuseMap model(m_map_file, 3.0);
        test_value(model.value(), 3.0);
        test_assert(model.filename() == m_map_file, "Expected \""+m_map_file+"\"");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test skymap value constructor
    test_try("Test skymap value constructor");
    try {
        GSkymap map("GAL", 16, "RING", 10);
        GModelSpatialDiffuseMap model(map, 3.0);
        test_value(model.value(), 3.0);
        test_assert(model.filename() == "", "Expected \"\"");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }
    
    // Test XML constructor and attribute methods
    test_try("Test XML constructor, value and gradients");
    try {
        // Test XML constructor
        GXml                    xml(m_xml_model_diffuse_map);
        GXmlElement*            element = xml.element(0)->element(0)->element("spatialModel", 0);
        GModelSpatialDiffuseMap model(*element);
        test_value(model.size(), 1);
        test_assert(model.type() == "SpatialMap", "Expected \"SpatialMap\"");
        test_value(model.value(), 1.0);
        test_assert(model.filename() == m_map_file,
                                        "Model filename \""+m_map_file+"\" expected.");

        // Test value method
        model.value(3.9);
        test_value(model.value(), 3.9);

        // Test load method
        model.load(m_map_file);
        test_assert(model.filename() == m_map_file,
                                        "Model filename \""+m_map_file+"\" expected.");

        // Test map method
        model.map(GSkymap("GAL", 16, "RING", 10));
        test_value(model.map().npix(), 3072);

        // Test operator access
        test_value(model["Prefactor"].value(), 3.9);
        test_value(model["Prefactor"].error(), 0.0);
        test_value(model["Prefactor"].gradient(), 0.0);
        model["Prefactor"].value(2.1);
        model["Prefactor"].error(1.9);
        model["Prefactor"].gradient(0.8);
        test_value(model["Prefactor"].value(), 2.1);
        test_value(model["Prefactor"].error(), 1.9);
        test_value(model["Prefactor"].gradient(), 0.8);

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
 * @brief Test GModelSpatialRadialDisk class
 ***************************************************************************/
void TestGModel::test_radial_disk(void)
{
    // Test void constructor
    test_try("Test void constructor");
    try {
        GModelSpatialRadialDisk model;
        test_assert(model.type() == "DiskFunction",
                                    "Model type \"DiskFunction\" expected.");
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
        test_assert(model.type() == "DiskFunction", "Expected \"DiskFunction\"");
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
        test_assert(model.type() == "GaussFunction",
                                    "Model type \"GaussFunction\" expected.");
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
        test_assert(model.type() == "GaussFunction", "Expected \"GaussFunction\"");
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
        test_assert(model.type() == "ShellFunction",
                                    "Model type \"ShellFunction\" expected.");
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
        test_assert(model.type() == "ShellFunction", "Expected \"ShellFunction\"");
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
 * @brief Test GModelSpectralConst class
 ***************************************************************************/
void TestGModel::test_const(void)
{
    // Test void constructor
    test_try("Test void constructor");
    try {
        GModelSpectralConst model;
        test_assert(model.type() == "ConstantValue",
                                    "Model type \"ConstantValue\" expected.");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test value constructor
    test_try("Test value constructor");
    try {
        GModelSpectralConst model(3.0);
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
        GXml                xml(m_xml_model_point_const);
        GXmlElement*        element = xml.element(0)->element(0)->element("spectrum", 0);
        GModelSpectralConst model(*element);
        test_value(model.size(), 1);
        test_assert(model.type() == "ConstantValue", "Expected \"ConstantValue\"");
        test_value(model.value(), 5.7e-16);

        // Test value method
        model.value(3.9e-16);
        test_value(model.value(), 3.9e-16);

        // Test operator access
        const char* strarray[] = {"Value"};
        for (int i = 0; i < 1; ++i) {
            std::string keyname(strarray[i]);
            model[keyname].remove_range(); // To allow setting of any value
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
 * @brief Test GModelSpectralPlaw class
 ***************************************************************************/
void TestGModel::test_plaw(void)
{
    // Test void constructor
    test_try("Test void constructor");
    try {
        GModelSpectralPlaw model;
        test_assert(model.type() == "PowerLaw",
                                    "Model type \"PowerLaw\" expected.");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test value constructor
    test_try("Test value constructor");
    try {
        GModelSpectralPlaw model(2.0, -2.1, GEnergy(100.0, "MeV"));
        test_value(model.prefactor(), 2.0);
        test_value(model.index(), -2.1);
        test_value(model.pivot().MeV(), 100.0);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }
    
    // Test XML constructor and value
    test_try("Test XML constructor, value and gradients");
    try {
        // Test XML constructor
        GXml               xml(m_xml_model_point_plaw);
        GXmlElement*       element = xml.element(0)->element(0)->element("spectrum", 0);
        GModelSpectralPlaw model(*element);
        test_value(model.size(), 3);
        test_assert(model.type() == "PowerLaw", "Expected \"PowerLaw\"");
        test_value(model.prefactor(), 5.7e-16);
        test_value(model.index(), -2.48);
        test_value(model.pivot().TeV(), 0.3);

        // Test prefactor method
        model["Prefactor"].remove_range(); // To allow setting of any value
        model.prefactor(3.9);
        test_value(model.prefactor(), 3.9);

        // Test index method
        model["Index"].remove_range(); // To allow setting of any value
        model.index(2.1);
        test_value(model.index(), 2.1);

        // Test pivot method
        model["PivotEnergy"].remove_range(); // To allow setting of any value
        model.pivot(GEnergy(10.0, "MeV"));
        test_value(model.pivot().MeV(), 10.0);

        // Test operator access
        const char* strarray[] = {"Prefactor", "Index", "PivotEnergy"};
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
 * @brief Test GModelSpectralPlaw2 class
 ***************************************************************************/
void TestGModel::test_plaw2(void)
{
    // Test void constructor
    test_try("Test void constructor");
    try {
        GModelSpectralPlaw2 model;
        test_assert(model.type() == "PowerLaw2",
                                    "Model type \"PowerLaw2\" expected.");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test value constructor
    test_try("Test value constructor");
    try {
        GModelSpectralPlaw2 model(2.0, -2.1, GEnergy(10.0, "MeV"), GEnergy(100.0, "MeV"));
        test_value(model.integral(), 2.0);
        test_value(model.index(), -2.1);
        test_value(model.emin().MeV(), 10.0);
        test_value(model.emax().MeV(), 100.0);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }
    
    // Test XML constructor and value
    test_try("Test XML constructor, value and gradients");
    try {
        // Test XML constructor
        GXml                xml(m_xml_model_point_plaw2);
        GXmlElement*        element = xml.element(0)->element(0)->element("spectrum", 0);
        GModelSpectralPlaw2 model(*element);
        test_value(model.size(), 4);
        test_assert(model.type() == "PowerLaw2", "Expected \"PowerLaw2\"");
        test_value(model.integral(), 1.0e-7);
        test_value(model.index(), -2.0);
        test_value(model.emin().MeV(), 100.0);
        test_value(model.emax().MeV(), 500000.0);

        // Test integral method
        model.integral(2.1e-7);
        test_value(model.integral(), 2.1e-7);

        // Test index method
        model.index(-2.3);
        test_value(model.index(), -2.3);

        // Test emin method
        model.emin(GEnergy(10.0, "MeV"));
        test_value(model.emin().MeV(), 10.0);

        // Test emax method
        model.emax(GEnergy(10.0, "MeV"));
        test_value(model.emax().MeV(), 10.0);

        // Test operator access
        const char* strarray[] = {"Integral", "Index", "LowerLimit", "UpperLimit"};
        for (int i = 0; i < 4; ++i) {
            std::string keyname(strarray[i]);
            model[keyname].remove_range(); // To allow setting of any value
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
 * @brief Test GModelSpectralExpPlaw class
 ***************************************************************************/
void TestGModel::test_eplaw(void)
{
    // Test void constructor
    test_try("Test void constructor");
    try {
        GModelSpectralExpPlaw model;
        test_assert(model.type() == "ExpCutoff",
                                    "Model type \"ExpCutoff\" expected.");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test value constructor
    test_try("Test value constructor");
    try {
        GModelSpectralExpPlaw model(2.0, -2.1, GEnergy(100.0, "MeV"), GEnergy(1.0, "GeV"));
        test_value(model.prefactor(), 2.0);
        test_value(model.index(), -2.1);
        test_value(model.pivot().MeV(), 100.0);
        test_value(model.cutoff().GeV(), 1.0);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }
    
    // Test XML constructor and value
    test_try("Test XML constructor, value and gradients");
    try {
        // Test XML constructor
        GXml                  xml(m_xml_model_point_eplaw);
        GXmlElement*          element = xml.element(0)->element(0)->element("spectrum", 0);
        GModelSpectralExpPlaw model(*element);
        test_value(model.size(), 4);
        test_assert(model.type() == "ExpCutoff", "Expected \"ExpCutoff\"");
        test_value(model.prefactor(), 5.7e-16);
        test_value(model.index(), -2.48);
        test_value(model.pivot().TeV(), 0.3);
        test_value(model.cutoff().TeV(), 1.0);

        // Test prefactor method
        model.prefactor(2.3e-16);
        test_value(model.prefactor(), 2.3e-16);

        // Test index method
        model.index(-2.6);
        test_value(model.index(), -2.6);

        // Test pivot method
        model.pivot(GEnergy(0.5, "TeV"));
        test_value(model.pivot().TeV(), 0.5);

        // Test cutoff method
        model.cutoff(GEnergy(10.0, "TeV"));
        test_value(model.cutoff().TeV(), 10.0);

        // Test operator access
        const char* strarray[] = {"Prefactor", "Index", "PivotEnergy", "Cutoff"};
        for (int i = 0; i < 4; ++i) {
            std::string keyname(strarray[i]);
            model[keyname].remove_range(); // To allow setting of any value
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
 * @brief Test GModelSpectralBrokenPlaw class
 ***************************************************************************/
void TestGModel::test_bplaw(void)
{
    // Test void constructor
    test_try("Test void constructor");
    try {
        GModelSpectralBrokenPlaw model;
        test_assert(model.type() == "BrokenPowerLaw",
                                    "Model type \"BrokenPowerLaw\" expected.");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test value constructor
    test_try("Test value constructor");
    try {
        GModelSpectralBrokenPlaw model(2.0, -2.1, GEnergy(100.0, "MeV"), -2.8);
        test_value(model.prefactor(), 2.0);
        test_value(model.index1(), -2.1);
        test_value(model.breakenergy().MeV(), 100.0);
        test_value(model.index2(), -2.8);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }
    
    // Test XML constructor and value
    test_try("Test XML constructor, value and gradients");
    try {
        // Test XML constructor
        GXml                     xml(m_xml_model_point_bplaw);
        GXmlElement*             element = xml.element(0)->element(0)->element("spectrum", 0);
        GModelSpectralBrokenPlaw model(*element);
        test_value(model.size(), 4);
        test_assert(model.type() == "BrokenPowerLaw", "Expected \"BrokenPowerLaw\"");
        test_value(model.prefactor(), 5.7e-16);
        test_value(model.index1(), -2.48);
        test_value(model.breakenergy().TeV(), 0.3);
        test_value(model.index2(), -2.70);

        // Test prefactor method
        model.prefactor(2.3e-16);
        test_value(model.prefactor(), 2.3e-16);

        // Test index1 method
        model.index1(-2.6);
        test_value(model.index1(), -2.6);

        // Test breakenergy method
        model.breakenergy(GEnergy(0.5, "TeV"));
        test_value(model.breakenergy().TeV(), 0.5);

        // Test index2 method
        model.index2(-3.6);
        test_value(model.index2(), -3.6);

        // Test operator access
        const char* strarray[] = {"Prefactor", "Index1", "BreakValue", "Index2"};
        for (int i = 0; i < 4; ++i) {
            std::string keyname(strarray[i]);
            model[keyname].remove_range(); // To allow setting of any value
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
 * @brief Test GModelSpectralLogParabola class
 ***************************************************************************/
void TestGModel::test_logparabola(void)
{
    // Test void constructor
    test_try("Test void constructor");
    try {
        GModelSpectralLogParabola model;
        test_assert(model.type() == "LogParabola",
                                    "Model type \"LogParabola\" expected.");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test value constructor
    test_try("Test value constructor");
    try {
        GModelSpectralLogParabola model(2.0, -2.1, GEnergy(100.0, "MeV"), -0.2);
        test_value(model.prefactor(), 2.0);
        test_value(model.index(), -2.1);
        test_value(model.pivot().MeV(), 100.0);
        test_value(model.curvature(), -0.2);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }
    
    // Test XML constructor and value
    test_try("Test XML constructor, value and gradients");
    try {
        // Test XML constructor
        GXml                      xml(m_xml_model_point_logparabola);
        GXmlElement*              element = xml.element(0)->element(0)->element("spectrum", 0);
        GModelSpectralLogParabola model(*element);
        test_value(model.size(), 4);
        test_assert(model.type() == "LogParabola", "Expected \"LogParabola\"");
        test_value(model.prefactor(), 5.878e-16);
        test_value(model.index(), -2.32473);
        test_value(model.pivot().TeV(), 1.0);
        test_value(model.curvature(), -0.074);

        // Test prefactor method
        model.prefactor(2.3e-16);
        test_value(model.prefactor(), 2.3e-16);

        // Test index method
        model.index(-2.6);
        test_value(model.index(), -2.6);

        // Test pivot method
        model.pivot(GEnergy(0.5, "TeV"));
        test_value(model.pivot().TeV(), 0.5);

        // Test curvature method
        model.curvature(-0.1);
        test_value(model.curvature(), -0.1);

        // Test operator access
        const char* strarray[] = {"Prefactor", "Index", "PivotEnergy", "Curvature"};
        for (int i = 0; i < 4; ++i) {
            std::string keyname(strarray[i]);
            model[keyname].remove_range(); // To allow setting of any value
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
 * @brief Test GModelSpectralNodes class
 ***************************************************************************/
void TestGModel::test_nodes(void)
{
    // Test void constructor
    test_try("Test void constructor");
    try {
        GModelSpectralNodes model;
        test_assert(model.type() == "NodeFunction",
                                    "Model type \"NodeFunction\" expected.");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test node function manimulation
    test_try("Test node function manimulation");
    try {
        GModelSpectralNodes model;
        model.reserve(3);
        test_value(model.size(), 0);
        test_value(model.nodes(), 0);
        model.append(GEnergy(1.0, "MeV"), 1.0);
        test_value(model.size(), 2);
        test_value(model.nodes(), 1);
        test_value(model.energy(0).MeV(), 1.0);
        test_value(model.intensity(0), 1.0);
        model.append(GEnergy(10.0, "MeV"), 0.1);
        test_value(model.size(), 4);
        test_value(model.nodes(), 2);
        test_value(model.energy(0).MeV(), 1.0);
        test_value(model.energy(1).MeV(), 10.0);
        test_value(model.intensity(0), 1.0);
        test_value(model.intensity(1), 0.1);
        model.remove(0);
        test_value(model.size(), 2);
        test_value(model.nodes(), 1);
        test_value(model.energy(0).MeV(), 10.0);
        test_value(model.intensity(0), 0.1);
        model.insert(0, GEnergy(1.0, "MeV"), 1.0);
        test_value(model.size(), 4);
        test_value(model.nodes(), 2);
        test_value(model.energy(0).MeV(), 1.0);
        test_value(model.energy(1).MeV(), 10.0);
        test_value(model.intensity(0), 1.0);
        test_value(model.intensity(1), 0.1);
        model.extend(model);
        test_value(model.size(), 8);
        test_value(model.nodes(), 4);
        test_value(model.energy(0).MeV(), 1.0);
        test_value(model.energy(1).MeV(), 10.0);
        test_value(model.energy(2).MeV(), 1.0);
        test_value(model.energy(3).MeV(), 10.0);
        test_value(model.intensity(0), 1.0);
        test_value(model.intensity(1), 0.1);
        test_value(model.intensity(2), 1.0);
        test_value(model.intensity(3), 0.1);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }
   
    // Test XML constructor and value
    test_try("Test XML constructor, value and gradients");
    try {
        // Test XML constructor
        GXml                      xml(m_xml_model_point_nodes);
        GXmlElement*              element = xml.element(0)->element(0)->element("spectrum", 0);
        GModelSpectralNodes model(*element);
        test_value(model.size(), 4);
        test_value(model.nodes(), 2);
        test_assert(model.type() == "NodeFunction", "Expected \"NodeFunction\"");
        test_value(model.energy(0).MeV(), 1.0);
        test_value(model.energy(1).MeV(), 10.0);
        test_value(model.intensity(0), 1.0e-7);
        test_value(model.intensity(1), 0.1e-7);

        // Test energy method
        model.energy(0, GEnergy(0.1, "MeV"));
        test_value(model.energy(0).MeV(), 0.1);

        // Test intensity method
        model.intensity(0, 2.0e-7);
        test_value(model.intensity(0), 2.0e-7);

        // Test operator access
        const char* strarray[] = {"Energy0", "Energy1", "Intensity0", "Intensity1"};
        for (int i = 0; i < 4; ++i) {
            std::string keyname(strarray[i]);
            model[keyname].remove_range(); // To allow setting of any value
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
 * @brief Test GModelSpectralFunc class
 ***************************************************************************/
void TestGModel::test_filefct(void)
{
    // Test void constructor
    test_try("Test void constructor");
    try {
        GModelSpectralFunc model;
        test_assert(model.type() == "FileFunction",
                                    "Model type \"FileFunction\" expected.");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test value constructor
    test_try("Test value constructor");
    try {
        GModelSpectralFunc model("data/filefunction.txt", 2.0);
        test_assert(model.filename() == "data/filefunction.txt", "Expected \"data/filefunction.txt\"");
        test_value(model.norm(), 2.0);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }
   
    // Test XML constructor and value
    test_try("Test XML constructor, value and gradients");
    try {
        // Test XML constructor
        GXml               xml(m_xml_model_point_filefct);
        GXmlElement*       element = xml.element(0)->element(0)->element("spectrum", 0);
        GModelSpectralFunc model(*element);
        test_value(model.size(), 1);
        test_assert(model.type() == "FileFunction", "Expected \"FileFunction\"");
        test_assert(model.filename() == "data/filefunction.txt", "Expected \"data/filefunction.txt\"");
        test_value(model.norm(), 1.0);

        // Test filename method
        model.filename("data/filefunction.txt");
        test_assert(model.filename() == "data/filefunction.txt", "Expected \"data/filefunction.txt\"");

        // Test norm method
        model.norm(3.0);
        test_value(model.norm(), 3.0);

        // Test operator access
        const char* strarray[] = {"Normalization"};
        for (int i = 0; i < 1; ++i) {
            std::string keyname(strarray[i]);
            model[keyname].remove_range(); // To allow setting of any value
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
 * @brief Test GModelTemporalConst class
 ***************************************************************************/
void TestGModel::test_temp_const(void)
{
    // Test void constructor
    test_try("Test void constructor");
    try {
        GModelTemporalConst model;
        test_assert(model.type() == "Constant",
                                    "Model type \"Constant\" expected.");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test value constructor
    test_try("Test value constructor");
    try {
        GModelTemporalConst model(3.0);
        test_value(model.norm(), 3.0);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }
    
    // Test value
    test_try("Test XML constructor, value and gradients");
    try {
        // Test XML constructor
        GModelTemporalConst model;

        // Test value method
        model.norm(3.9);
        test_value(model.norm(), 3.9);

        // Test operator access
        const char* strarray[] = {"Constant"};
        for (int i = 0; i < 1; ++i) {
            std::string keyname(strarray[i]);
            model[keyname].remove_range(); // To allow setting of any value
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
 * @brief Test XML model.
 *
 * @param[in] name Model name.
 * @param[in] filename XML model filename.
 *
 * Test loading and saving of XML model.
 ***************************************************************************/
void TestGModel::test_xml_model(const std::string& name,
                                const std::string& filename)
{
    // Test load constructor
    test_try("Test load constructor");
    try {
        GModels models(filename);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test saving and reloading
    test_try("Test saving and reloading");
    try {
        GModels models(filename);
        models.save("test.xml");
        models.load("test.xml");
        models.save("test.xml");
        models.load("test.xml");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Exit test
    return;

}


/***********************************************************************//**
 * @brief Test spatial model XML reading and writing
 ***************************************************************************/
void TestGModel::test_spatial_model(void)
{
    // Test spatial models XML interface
    test_xml_model("GModelSpatialPointSource",    m_xml_model_point_plaw);
    test_xml_model("GModelSpatialRadialDisk",     m_xml_model_radial_disk);
    test_xml_model("GModelSpatialRadialGauss",    m_xml_model_radial_gauss);
    test_xml_model("GModelSpatialRadialShell",    m_xml_model_radial_shell);
    test_xml_model("GModelSpatialEllipticalDisk", m_xml_model_elliptical_disk);
    test_xml_model("GModelSpatialDiffuseConst",   m_xml_model_diffuse_const);
    test_xml_model("GModelSpatialDiffuseMap",     m_xml_model_diffuse_map);
    test_xml_model("GModelSpatialDiffuseCube",    m_xml_model_diffuse_cube);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test spectral model XML reading and writing
 ***************************************************************************/
void TestGModel::test_spectral_model(void)
{
    // Test spectral models XML interface
    test_xml_model("GModelSpectralConst",       m_xml_model_point_const);
    test_xml_model("GModelSpectralPlaw",        m_xml_model_point_plaw);
    test_xml_model("GModelSpectralPlaw2",       m_xml_model_point_plaw2);
    test_xml_model("GModelSpectralExpPaw",      m_xml_model_point_eplaw);
    test_xml_model("GModelSpectralBrokenPlaw",  m_xml_model_point_bplaw);
    test_xml_model("GModelSpectralLogParabola", m_xml_model_point_logparabola);
    test_xml_model("GModelSpectralNodes",       m_xml_model_point_nodes);
    test_xml_model("GModelSpectralFunc",        m_xml_model_point_filefct);

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
    // Test GModelRegistry void constructor
    test_try("Test GModelRegistry void constructor");
    try {
        GModelRegistry registry;
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

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

    // Test GModelSpatialRegistry void constructor
    test_try("Test GModelSpatialRegistry void constructor");
    try {
        GModelSpatialRegistry registry;
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test GModelSpatialRegistry allocators
    test_try("Test GModelSpatialRegistry model allocation");
    try {
        // Loop over all models in registry
        GModelSpatialRegistry registry;
        int num = registry.size();
        for (int i = 0; i < num; ++i) {
            GModelSpatial* ptr = registry.alloc(registry.name(i));
            test_assert(ptr != NULL, "Model pointer for \""+ registry.name(i)+"\" is NULL");
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

    // Test GModelSpectralRegistry void constructor
    test_try("Test GModelSpectralRegistry void constructor");
    try {
        GModelSpectralRegistry registry;
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test GModelSpectralRegistry allocators
    test_try("Test GModelSpectralRegistry model allocation");
    try {
        // Loop over all models in registry
        GModelSpectralRegistry registry;
        int num = registry.size();
        for (int i = 0; i < num; ++i) {
            GModelSpectral* ptr = registry.alloc(registry.name(i));
            test_assert(ptr != NULL, "Model pointer for \""+registry.name(i)+"\" is NULL");
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

    // Test GModelTemporalRegistry void constructor
    test_try("Test GModelTemporalRegistry void constructor");
    try {
        GModelTemporalRegistry registry;
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test GModelTemporalRegistry allocators
    test_try("Test GModelTemporalRegistry model allocation");
    try {
        // Loop over all models in registry
        GModelTemporalRegistry registry;
        int num = registry.size();
        for (int i = 0; i < num; ++i) {
            GModelTemporal* ptr = registry.alloc(registry.name(i));
            test_assert(ptr != NULL, "Model pointer for \""+ \
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
    bool was_successful=true;

    // Create a test suite
    TestGModel test;

    // Append to the container
    testsuite.append(test);

    // Run
    was_successful=testsuite.run();

    // Save xml report
    testsuite.save("reports/GModel.xml");

    // Return
    return was_successful ? 0:1;
}
