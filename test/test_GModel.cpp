/***************************************************************************
 *                  test_GModel.cpp  -  test GModel class                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

/* __ Includes ___________________________________________________________ */
#include <stdlib.h>
#include "test_GModel.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Globals ____________________________________________________________ */


/***************************************************************************
 *  Test: Model parameter handling                                         *
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


/***************************************************************************
 *  Test: Model handling                                                   *
 ***************************************************************************/
void test_model(void)
{
    // Write header
    std::cout << "Test GModel: ";

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
    GModel crab;
    try {
        crab = GModel(point_source, power_law);
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
        models.add(crab);
        models.add(crab);
        models.add(crab);
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
    GData data;
    GLATObservation obs;
    try {
        GModels models;
        models.add(crab);
        obs.load_unbinned("data/lat/ft1.fits.gz", "data/lat/ft2.fits.gz", "");
        data.add(obs);
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

    // Plot final test success
    std::cout << " ok." << std::endl;

    // Exit test
    return;

}


/***************************************************************************
 *                            Main test function                           *
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

    // Return
    return 0;
}
