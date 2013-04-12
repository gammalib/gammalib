/***************************************************************************
 *              test_GOptimizer.cpp - test GOptimizer class                *
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
 * @file test_GOptimizer.cpp
 * @brief Implementation of unit tests for optimizer module
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "test_GOptimizer.hpp"
#include "testinst/GTestLib.hpp"

/* __ Coding definitions _________________________________________________ */
#define RATE      13.0        //!< Events per seconde. For events generation.
#define UN_BINNED 0
#define BINNED    1


/***********************************************************************//**
 * @brief Set parameters and tests
 **************************************************************************/
void TestGOptimizer::set(void){

    // Test name
    name("Optimizer module");

    // Append tests
    append(static_cast<pfunction>(&TestGOptimizer::test_unbinned_optimizer), "Test unbinned optimization");
    append(static_cast<pfunction>(&TestGOptimizer::test_binned_optimizer), "Test binned optimization");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test optimizer
 *
 * @param[in] mode Testing mode.
 * 
 * This method supports two testing modes: 0 = unbinned and 1 = binned.
 ***************************************************************************/
void TestGOptimizer::test_optimizer(const int& mode)
{
    // Create Test Model
    GTestModelData model;

    // Create Models conteners
    GModels models;
    models.append(model);

    // Set time iterval
    GTime tmin(0.0);
    GTime tmax(1800.0);

    // Rate : events/sec
    double rate = RATE;

    // Create observations
    GObservations obs;

    // Add some observation
    for (int i = 0; i < 6; ++i) {

        // Random Generator
        GRan ran;
        ran.seed(i);

        // Allocate events pointer
        GEvents *events;

        // Create either a event list or an event cube
        if (mode == UN_BINNED) {
            events = model.generateList(rate,tmin,tmax,ran);
        }
        else {
            events = model.generateCube(rate,tmin,tmax,ran);
        }

        // Create an observation
        GTestObservation ob;
        ob.id(gammalib::str(i));

        // Add events to the observation
        ob.events(events);
        ob.ontime(tmax.secs()-tmin.secs());
        obs.append(ob);
    }

    // Add the model to the observation
    obs.models(models);

    // Create a GLog for show the interations of optimizer.
    GLog log;

    // Create an optimizer.
    GOptimizerLM opt(log);

    opt.max_stalls(50);

    // Optimize
    obs.optimize(opt);

    // Get the result
    GModelPar result = (*(obs.models()[0]))[0];

    // Check if converged
    test_assert(opt.status()==0, "Check if converged", 
                                 "Optimizer did not converge"); 

    // Check if value is correct
    test_value(result.factor_value(), RATE, result.factor_error()*3); 

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test unbinned optimizer
 ***************************************************************************/
void TestGOptimizer::test_unbinned_optimizer(void)
{
    // Test
    test_optimizer(UN_BINNED);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test binned optimizer
 ***************************************************************************/
void TestGOptimizer::test_binned_optimizer(void)
{
    // Test
    test_optimizer(BINNED);

    // Return
    return;
}


/***************************************************************************
 * @brief Main entry point for test executable
 ***************************************************************************/
int main(void)
{
    // Allocate test suit container
    GTestSuites testsuites("Optimizer module");

    // Initially assume that we pass all tests
    bool success = true;

    // Create and append test suite
    TestGOptimizer openmp;
    testsuites.append(openmp);

    // Run the testsuites
    success = testsuites.run();

    // Save test report
    testsuites.save("reports/GOptimizer.xml");

    // Return success status
    return (success ? 0 : 1);
}
