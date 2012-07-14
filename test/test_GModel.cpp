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
void test_model_par(void)
{
    // Write header
    std::cout << "Test GModelPar: ";

    // Load unbinned LAT observation
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
    }
    catch (std::exception &e) {
        std::cout << std::endl << "TEST ERROR: Unable to handle model parameter." 
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Plot final test success
    std::cout << " ok." << std::endl;

    // Exit test
    return;

}


/***********************************************************************//**
 * @brief Test model handling.
 ***************************************************************************/
void test_model(void)
{
    // Write header
    std::cout << "Test GModelPointSource: ";

    // Setup spatial model
    GModelSpatialPtsrc point_source;
    try {
        GSkyDir dir;
        dir.radec_deg(83.6331, +22.0145);
        point_source = GModelSpatialPtsrc(dir);
    }
    catch (std::exception &e) {
        std::cout << std::endl << "TEST ERROR: Unable setup GModelSpatialPtsrc." 
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";
    if (point_source.ra() != 83.6331 || point_source.dec() != +22.0145) {
        std::cout << std::endl << "TEST ERROR: Bad values in GModelSpatialPtsrc."
                  << std::endl;
        throw;
    }
    std::cout << ".";
    //std::cout << point_source << std::endl;

    // Setup spectral model
    GModelSpectralPlaw power_law;
    try {
        power_law = GModelSpectralPlaw(1.0e-7, -2.1);
    }
    catch (std::exception &e) {
        std::cout << std::endl << "TEST ERROR: Unable setup GModelSpectralPlaw." 
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";
    if (power_law.norm() != 1.0e-7 || power_law.index() != -2.1) {
        std::cout << std::endl << "TEST ERROR: Bad values in GModelSpectralPlaw."
                  << std::endl;
        throw;
    }
    std::cout << ".";
    //std::cout << power_law << std::endl;

    // Setup Crab model
    GModelPointSource crab;
    try {
        crab = GModelPointSource(point_source, power_law);
        crab.name("Crab");
    }
    catch (std::exception &e) {
        std::cout << std::endl << "TEST ERROR: Unable setup GModel." 
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";
    //std::cout << crab << std::endl;

    // Put model in container
    GModels models;
    try {
        models.append(crab);
        models.append(crab);
        models.append(crab);
    }
    catch (std::exception &e) {
        std::cout << std::endl << "TEST ERROR: Unable setup GModels." 
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";
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
    // Plot final test success
    std::cout << " ok." << std::endl;

    // Exit test
    return;

}


/***********************************************************************//**
 * @brief Test XML model.
 ***************************************************************************/
void test_xml_model(const std::string &name, const std::string &filename)
{
    // Write header
    std::cout << "Test " << name << ": ";

    // Test load constructor
    try {
        GModels models(filename);
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to load model from XML document."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Test save constructor
    try {
        GModels models(filename);
        models.save("test.xml");
        models.load("test.xml");
        models.save("test.xml");
        models.load("test.xml");
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to save and load model from XML document."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Plot final test success
    std::cout << " ok." << std::endl;

    // Exit test
    return;

}


/***********************************************************************//**
 * @brief Test models.
 ***************************************************************************/
void test_models(void)
{
    // Write header
    std::cout << "Test GModels: ";

    // Test void constructor
    try {
        GModels models;
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to construct empty model container."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Test load constructor
    try {
        GModels models(xml_file);
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to construct model container from XML document."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Test saving and loading
    try {
        GModels models(xml_file);
        models.save("test.xml");
        models.load("test.xml");
        models.save("test.xml");
        models.load("test.xml");
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to save model container in XML document."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Plot final test success
    std::cout << " ok." << std::endl;

    // Exit test
    return;

}


/***********************************************************************//**
 * @brief Main test function.
 ***************************************************************************/
int main(void)
{
    // Dump header
    std::cout << std::endl;
    std::cout << "************************" << std::endl;
    std::cout << "* GModel class testing *" << std::endl;
    std::cout << "************************" << std::endl;

    // Execute the tests
    test_model_par();
    test_model();
    test_models();

    // Test spectral models
    test_xml_model("GModelSpectralNodes", xml_model_point_nodes);
    test_xml_model("GModelSpatialMap", xml_model_spatial_map);

    // Return
    return 0;
}
