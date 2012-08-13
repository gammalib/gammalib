/***************************************************************************
 *                  test_GModel.cpp  -  test GModel class                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2012 by Juergen Knoedlseder                         *
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
 * @author J. Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iostream>
#include <ostream>
#include <stdexcept>
#include <stdlib.h>
#include "GammaLib.hpp"

/* __ Globals ____________________________________________________________ */
const std::string xml_file              = "data/crab.xml";
const std::string xml_model_point_nodes = "data/model_point_nodes.xml";
const std::string xml_model_spatial_map = "data/model_spatial_map.xml";


/***********************************************************************//**
 * @brief Test model parameter handling.
 ***************************************************************************/
void test_model_par(GTestSuite& testsuite)
{
    // Load unbinned LAT observation
    testsuite.test_try("Load unbinned LAT observation");
    try {
        GModelPar par;
        par.value(47.01);
        par.error(2.003);
        par.name("Test parameter");
        par.unit("MeV");
//        std::cout << par << std::endl;
        par.fix();
//        std::cout << par << std::endl;
        par.min(2.0);
//        std::cout << par << std::endl;
        par.max(200.0);
//        std::cout << par << std::endl;
        par.remove_min();
//        std::cout << par << std::endl;
        par.remove_max();
//        std::cout << par << std::endl;

        testsuite.test_try_success();
    }
    catch (std::exception &e) {
        testsuite.test_try_failure(e);
    }

    // Exit test
    return;

}


/***********************************************************************//**
 * @brief Test model handling.
 ***************************************************************************/
void test_model(GTestSuite& testsuite)
{
    // Setup spatial model
    GModelSpatialPtsrc point_source;
    testsuite.test_try("Setup spatial model");
    try {
        GSkyDir dir;
        dir.radec_deg(83.6331, +22.0145);
        point_source = GModelSpatialPtsrc(dir);

        testsuite.test_try_success();
    }
    catch (std::exception &e) {
        testsuite.test_try_failure(e);
    }

    testsuite.test_assert((point_source.ra() == 83.6331 && point_source.dec() == +22.0145),
                           "Test if ra=83.6331 and dec=22.0145",
                           "Bad values in GModelSpatialPtsrc");

    //std::cout << point_source << std::endl;

    // Setup spectral model
    GModelSpectralPlaw power_law;
    testsuite.test_try("Setup spectral model");
    try {
        power_law = GModelSpectralPlaw(1.0e-7, -2.1);
        testsuite.test_try_success();
    }
    catch (std::exception &e) {
        testsuite.test_try_failure(e);
    }

    testsuite.test_assert((power_law.norm() == 1.0e-7 && power_law.index() == -2.1),
                           "Test if norm=1.0e-7 and index=-2.1","Bad values in GModelSpectralPlaw");

    //std::cout << power_law << std::endl;

    // Setup Crab model
    GModelPointSource crab;
    testsuite.test_try("Setup Crab model");
    try {
        crab = GModelPointSource(point_source, power_law);
        crab.name("Crab");

        testsuite.test_try_success();
    }
    catch (std::exception &e) {
        testsuite.test_try_failure(e);
    }

    //std::cout << crab << std::endl;

    // Put model in container
    GModels models;
    testsuite.test_try("Put model in container");
    try {
        models.append(crab);
        models.append(crab);
        models.append(crab);

        testsuite.test_try_success();
    }
    catch (std::exception &e) {
        testsuite.test_try_failure(e);
    }

    //std::cout << models << std::endl;

    // Put container in data
/*
    GObservations   data;
    GLATObservation obs;
    try {
        GModels models;
        models.append(crab);
        obs.load_unbinned("data/lat/ft1.fits.gz", "data/lat/ft2.fits.gz", "");
        data.append(obs);
        data.models(models);
    }
    catch (std::exception &e) {
        std::cout << std::endl << "TEST ERROR: Unable setup GData with GModels." 
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";
    //std::cout << data << std::endl;
*/

    // Exit test
    return;

}


/***********************************************************************//**
 * @brief Test XML model.
 ***************************************************************************/
void test_xml_model(const std::string &name, const std::string &filename,GTestSuite& testsuite)
{
    // Test load constructor
    testsuite.test_try("Test load constructor");
    try {
        GModels models(filename);
        testsuite.test_try_success();
    }
    catch (std::exception &e) {
        testsuite.test_try_failure(e);
    }

    // Test save constructor
    testsuite.test_try("Test load constructor");
    try {
        GModels models(filename);
        models.save("test.xml");
        models.load("test.xml");
        models.save("test.xml");
        models.load("test.xml");

        testsuite.test_try_success();
    }
    catch (std::exception &e) {
        testsuite.test_try_failure(e);
    }

    // Exit test
    return;

}


/***********************************************************************//**
 * @brief Test models.
 ***************************************************************************/
void test_models(GTestSuite& testsuite)
{
    // Test void constructor
    testsuite.test_try("Test void constructor");
    try {
        GModels models;
        testsuite.test_try_success();
    }
    catch (std::exception &e) {
        testsuite.test_try_failure(e);
    }

    // Test load constructor
    testsuite.test_try("Test load constructor");
    try {
        GModels models(xml_file);
        testsuite.test_try_success();
    }
    catch (std::exception &e) {
        testsuite.test_try_failure(e);
    }

    // Test saving and loading
    testsuite.test_try("Test saving and loading");
    try {
        GModels models(xml_file);
        models.save("test.xml");
        models.load("test.xml");
        models.save("test.xml");
        models.load("test.xml");

        testsuite.test_try_success();
    }
    catch (std::exception &e) {
        testsuite.test_try_failure(e);
    }

    // Exit test
    return;

}


/***********************************************************************//**
 * @brief Main test function.
 ***************************************************************************/
int main(void)
{
    //Create a test suites
    GTestSuites testsuites("GModel class testing");

    // Create a test suite
    GTestSuite * testsuite = NULL;

    // Test GModelPar
    testsuite = new GTestSuite("GModelPar");
    testsuites.append(*testsuite);
    test_model_par(*testsuite);
    testsuite->end_test();

    // Test GModel
    testsuite = new GTestSuite("GModelPointSource");
    testsuites.append(*testsuite);
    test_model(*testsuite);
    testsuite->end_test();

    // Test GModels
    testsuite = new GTestSuite("GModels");
    testsuites.append(*testsuite);
    test_models(*testsuite);
    testsuite->end_test();

    // Test spectral models

    //Test GModelSpectraNodes
    testsuite = new GTestSuite("GModelSpectralNodes");
    testsuites.append(*testsuite);
    test_xml_model("GModelSpectralNodes", xml_model_point_nodes,*testsuite);
    testsuite->end_test();

    //Test GModelSpatialMap
    testsuite = new GTestSuite("GModelSpatialMap");
    testsuites.append(*testsuite);
    test_xml_model("GModelSpatialMap", xml_model_spatial_map,*testsuite);
    testsuite->end_test();

    //save xml report
    testsuites.save("reports/GModel.xml");

    // Return
    return 0;
}
