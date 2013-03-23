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
    m_xml_file                  = "data/crab.xml";
    m_xml_model_point_nodes     = "data/model_point_nodes.xml";
    m_xml_model_diffuse_const   = "data/model_spatial_const.xml";
    m_xml_model_diffuse_map     = "data/model_spatial_map.xml";
    m_xml_model_radial_disk     = "data/model_radial_disk.xml";
    m_xml_model_radial_gauss    = "data/model_radial_gauss.xml";
    m_xml_model_radial_shell    = "data/model_radial_shell.xml";
    m_xml_model_elliptical_disk = "data/model_elliptical_disk.xml";

    // Add tests
    add_test(static_cast<pfunction>(&TestGModel::test_model_par), "Test GModelPar");
    add_test(static_cast<pfunction>(&TestGModel::test_sky_model), "Test GModelSky");
    add_test(static_cast<pfunction>(&TestGModel::test_diffuse_const), "Test GModelSpatialDiffuseConst");
    add_test(static_cast<pfunction>(&TestGModel::test_spatial_model), "Test spatial model XML I/O");
    add_test(static_cast<pfunction>(&TestGModel::test_spectral_model), "Test spectral model XML I/O");
    //
    add_test(static_cast<pfunction>(&TestGModel::test_model), "Test model handling");
    add_test(static_cast<pfunction>(&TestGModel::test_models), "Test models");

    // Return
    return;
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
        test_assert(!par.hasmin(), "Parameter shall have no minimum.");
        test_assert(!par.hasmax(), "Parameter shall have no maximum.");
        test_assert(!par.hasrange(), "Parameter shall have no range.");
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
        test_assert(par.isfree(), "Parameter freezing");
        test_assert(par.print() == 
                    "  Test parameter ...........: 94.02 +/- 4.006 MeV"
                    " (free,scale=2)", "Parameter printing");
        GModelPar par2("Another test parameter", 3.14, 3.0);
        test_value(par2.value(), 9.42);
        test_value(par2.factor_value(), 3.14);
        test_value(par2.scale(), 3.0);
        test_assert(par2.name() == "Another test parameter", "Parameter name");
        test_assert(par2.isfree(), "Parameter freezing");
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
        test_assert(!par.hasmin(), "Parameter shall have no minimum.");
        test_assert(!par.hasmax(), "Parameter shall have no maximum.");
        test_assert(!par.hasrange(), "Parameter shall have no range.");
        par.min(0.5);
        test_assert(par.hasmin(), "Parameter shall have minimum.");
        test_assert(!par.hasmax(), "Parameter shall have no maximum.");
        test_assert(!par.hasrange(), "Parameter shall have no range.");
        par.max(2.0);
        test_assert(par.hasmin(), "Parameter shall have minimum.");
        test_assert(par.hasmax(), "Parameter shall have maximum.");
        test_assert(par.hasrange(), "Parameter shall have range.");
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
        test_value(sky.value(dir, energy, time), 1.73e-07);

        // Test gradient
        GVector vector = sky.gradients(dir, energy, time);
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
 * @brief Test model handling.
 ***************************************************************************/
void TestGModel::test_model(void)
{
    // Setup spatial model
    GModelSpatialPointSource point_source;
    test_try("Setup spatial model");
    try {
        GSkyDir dir;
        dir.radec_deg(83.6331, +22.0145);
        point_source = GModelSpatialPointSource(dir);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    test_assert((point_source.ra() == 83.6331 && point_source.dec() == +22.0145),
                "Test if ra=83.6331 and dec=22.0145",
                "Bad values in GModelSpatialPointSource");

    // Setup spectral model
    GModelSpectralPlaw power_law;
    test_try("Setup spectral model");
    try {
        power_law = GModelSpectralPlaw(1.0e-7, -2.1, 100.0);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    test_assert((power_law.norm() == 1.0e-7 && power_law.index() == -2.1),
                "Test if norm=1.0e-7 and index=-2.1",
                "Bad values in GModelSpectralPlaw");

    // Setup Crab model
    GModelSky crab;
    test_try("Setup Crab model");
    try {
        crab = GModelSky(point_source, power_law);
        crab.name("Crab");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Put model in container
    GModels models;
    test_try("Put model in container");
    try {
        models.append(crab);
        models.append(crab);
        models.append(crab);
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
    GModels models2;
    models2.append(crab);
    models2.save("test_instrument.xml");
    models2.load("test_instrument.xml");
    test_value(models2[0]->scale("LAT").value(), 1.0);
    test_value(models2[0]->scale("CTA").value(), 0.5);
    test_value(models2[0]->scale("COM").value(), 1.0);

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
    test_xml_model("GModelSpatialDiffuseConst",   m_xml_model_diffuse_const);
    test_xml_model("GModelSpatialDiffuseMap",     m_xml_model_diffuse_map);
    test_xml_model("GModelSpatialRadialDisk",     m_xml_model_radial_disk);
    test_xml_model("GModelSpatialRadialGauss",    m_xml_model_radial_gauss);
    test_xml_model("GModelSpatialRadialShell",    m_xml_model_radial_shell);
    test_xml_model("GModelSpatialEllipticalDisk", m_xml_model_elliptical_disk);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test spectral model XML reading and writing
 ***************************************************************************/
void TestGModel::test_spectral_model(void)
{
    // Test spectral models XML interface
    test_xml_model("GModelSpectralNodes", m_xml_model_point_nodes);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test models.
 ***************************************************************************/
void TestGModel::test_models(void)
{
    // Test void constructor
    test_try("Test void constructor");
    try {
        GModels models;
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

    // Exit test
    return;

}


/***********************************************************************//**
 * @brief Main test function.
 ***************************************************************************/
int main(void)
{
    GTestSuites testsuite("GModel");

    bool was_successful=true;

    //Create a test suite
    TestGModel test;

    //Append to the container
    testsuite.append(test);

    //Run
    was_successful=testsuite.run();

    //save xml report
    testsuite.save("reports/GModel.xml");

    // Return
    return was_successful ? 0:1;
}
