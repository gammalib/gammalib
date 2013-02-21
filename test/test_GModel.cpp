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
    m_xml_file               = "data/crab.xml";
    m_xml_model_point_nodes  = "data/model_point_nodes.xml";
    m_xml_model_spatial_map  = "data/model_spatial_map.xml";
    m_xml_model_radial_disk  = "data/model_radial_disk.xml";
    m_xml_model_radial_gauss = "data/model_radial_gauss.xml";
    m_xml_model_radial_shell = "data/model_radial_shell.xml";

    // Add tests
    add_test(static_cast<pfunction>(&TestGModel::test_model_par), "Test model parameter handling");
    add_test(static_cast<pfunction>(&TestGModel::test_model), "Test model handling");
    add_test(static_cast<pfunction>(&TestGModel::test_models), "Test models");
    add_test(static_cast<pfunction>(&TestGModel::test_spectral_model), "Test spectral model");
    add_test(static_cast<pfunction>(&TestGModel::test_spatial_model), "Test spatial model");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test model parameter handling.
 *
 * Test setting and reading back of parameter value attributes.
 ***************************************************************************/
void TestGModel::test_model_par(void)
{
    // Set model parameter
    GModelPar par;
    par.value(47.01);
    par.error(2.003);
    par.scale(2.0);
    par.name("Test parameter");
    par.unit("MeV");
    par.free();
    par.min(2.0);
    par.max(200.0);
    par.remove_min();
    par.remove_max();

    // Check parameter access
    test_value(par.real_value(), 94.02);
    test_value(par.value(), 47.01);
    test_value(par.real_error(), 4.006);
    test_value(par.error(), 2.003);
    test_value(par.scale(), 2.0);
    test_assert(par.name() == "Test parameter", "Parameter name");
    test_assert(par.unit() == "MeV", "Parameter unit");
    test_assert(par.isfree(), "Parameter freezing");
    test_assert(par.print() == 
                "  Test parameter ...........: 94.02 +/- 4.006 MeV (free,scale=2)",
                "Parameter printing");

    // Set model parameter
    GModelPar par2("Another test parameter", 3.14, 3.0);
    test_value(par2.real_value(), 9.42);
    test_value(par2.value(), 3.14);
    test_value(par2.scale(), 3.0);
    test_assert(par2.name() == "Another test parameter", "Parameter name");
    test_assert(par2.isfree(), "Parameter freezing");
    test_assert(par2.print() == 
                "  Another test parameter ...: 9.42 +/- 0  (free,scale=3)",
                "Parameter printing");

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
        power_law = GModelSpectralPlaw(1.0e-7, -2.1);
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
    test_value(crab.scale("LAT").real_value(), 1.0);
    test_value(crab.scale("CTA").real_value(), 0.5);
    test_value(crab.scale("COM").real_value(), 1.0);

    // Test saving and loading
    GModels models2;
    models2.append(crab);
    models2.save("test_instrument.xml");
    models2.load("test_instrument.xml");
    test_value(models2[0].scale("LAT").real_value(), 1.0);
    test_value(models2[0].scale("CTA").real_value(), 0.5);
    test_value(models2[0].scale("COM").real_value(), 1.0);

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
 * @brief Test spectral model
 ***************************************************************************/
void TestGModel::test_spectral_model(void)
{
    // Test spectral models XML interface
    test_xml_model("GModelSpectralNodes", m_xml_model_point_nodes);

    // Return
    return;
}

/***********************************************************************//**
 * @brief Test spacial model
 ***************************************************************************/
void TestGModel::test_spatial_model(void)
{
    // Test spatial models XML interface
    test_xml_model("GModelSpatialMap",  m_xml_model_spatial_map);
    test_xml_model("GModelRadialDisk",  m_xml_model_radial_disk);
    test_xml_model("GModelRadialGauss", m_xml_model_radial_gauss);
    test_xml_model("GModelRadialShell", m_xml_model_radial_shell);

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
