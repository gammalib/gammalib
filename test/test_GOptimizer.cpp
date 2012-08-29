/***************************************************************************
 *             test_GOptimizer.cpp  -  test GOptimizer class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009 by Jurgen Knodlseder                                *
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

/* __ Includes ___________________________________________________________ */
#include <stdlib.h>
#include "test_GOptimizer.hpp"

/***************************************************************************
 *  Set                                                                    *
 ***************************************************************************/
void TestGOptimizer::set(void){
    // Test name
    name("GOptimizer");

    //add tests
    add_test(static_cast<pfunction>(&TestGOptimizer::test_optimizer),"Test LAT Response");

    return;
}

/***************************************************************************
 *  Test: Optimizer                                                        *
 ***************************************************************************/
void TestGOptimizer::test_optimizer(void)
{

    // Number of observations in data
    int nobs = 1;

    // Setup GData for optimizing
    GData           data;
    GLATObservation obs;
    test_try("Setup GData for optimizing");
    try {
        obs.load_binned("data/lat/cntmap.fits.gz", "", "");
        for (int i = 0; i < nobs; ++i)
            data.add(obs);

        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    //std::cout << data << std::endl;

    // Setup GModels for optimizing
    GModelSpatialPtsrc point_source;
    GModelSpectralPlaw power_law;
    GModel             crab;
    GModels            models;
    test_try("Setup GModels for optimizing");
    try {
        GSkyDir dir;
        dir.radec_deg(83.6331, +22.0145);
        point_source = GModelSpatialPtsrc(dir);
        power_law    = GModelSpectralPlaw(1.0e-7, -2.1);
        crab         = GModel(point_source, power_law);
        crab.name("Crab");
        models.add(crab);
        data.models(models);

        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    //std::cout << models << std::endl;

    // Setup parameters for optimizing
    test_try("Setup parameters for optimizing");
    GOptimizerPars pars;
    try {
        pars = GOptimizerPars(models);

        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    //std::cout << pars << std::endl;

    // Time simple GData iterator
    double  t_elapse;
    test_try("Time simple GData iterator");
    try {
        clock_t t_start = clock();
        int num = 0;
        int sum = 0;
        for (GData::iterator event = data.begin(); event != data.end(); ++event) {
            num++;
            sum += (int)event->counts();
        }
        t_elapse = double(clock() - t_start)/double(CLOCKS_PER_SEC);

        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    //std::cout << " - Reference time for GData::iterator: " << t_elapse << std::endl;

    // Setup LM optimizer
    GOptimizerLM opt;
    test_try("Setup LM optimizer");
    try {
        data.optimize(opt);
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    //std::cout << opt << std::endl;
    //std::cout << data << std::endl;

    // Exit test
    return;

}


/***************************************************************************
 *                            Main test function                           *
 ***************************************************************************/
int main(void)
{
    GTestSuites testsuites("GOptimizer");

    bool was_successful=true;

    //Create a test suite
    TestGOptimizer test;

    //Append to the container
    testsuites.append(test);

    //Run
    was_successful=testsuites.run();

    //save xml report
    testsuites.save("reports/GOptimizer.xml");

    // Return
    return was_successful ? 0:1;
}
