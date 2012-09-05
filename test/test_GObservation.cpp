/***************************************************************************
 *             test_GObservation.cpp  -  Test observation module           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012 by Jean-Baptiste Cayrou                             *
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
 * @file test_GObservation.cpp
 * @brief Test observation module
 * @author J.-B. Cayrou
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iostream>
#include "testinst/GTestLib.hpp"
#include "test_GObservation.hpp"

/* __ OpenMP section _____________________________________________________ */
#ifdef _OPENMP
#include <omp.h>
#ifdef __APPLE__ 
#ifdef __MACH__
#include <pthread.h>
pthread_attr_t gomp_thread_attr;
#endif
#endif
#endif

/* __ Coding definitions _________________________________________________ */
#define RATE      13.0        //!< Events per seconde. For events generation.
#define UN_BINNED 0
#define BINNED    1


#ifdef _OPENMP
/***********************************************************************//**
* @brief Set tests
***************************************************************************/
void TestOpenMP::set(void)
{
    // Set test name
    name("OpenMP");

    // Unbinned
    append(static_cast<pfunction>(&TestOpenMP::test_observations_optimizer_unbinned_1), "Test unbinned optimization (1 thread)");
    append(static_cast<pfunction>(&TestOpenMP::test_observations_optimizer_unbinned_10), "Test unbinned optimization (10 threads)");

    // Binned
    append(static_cast<pfunction>(&TestOpenMP::test_observations_optimizer_binned_1), "Test binned optimization (1 thread)");
    append(static_cast<pfunction>(&TestOpenMP::test_observations_optimizer_binned_10), "Test binned optimisation (10 threads)");

    // Return
    return;
}


/***********************************************************************//**
* @brief Test observations optimizer.
*
* @param[in] mode Testing mode.
* 
* This method supports two testing modes: 0 = unbinned and 1 = binned.
***************************************************************************/
GModelPar& TestOpenMP::test_observations_optimizer(int mode)
{
    // Create Test Model
    GTestModelData model;

    // Create Models conteners
    GModels models;
    models.append(model);

    // Time iterval
    GTime tmin(0,0,   "sec");
    GTime tmax(1800,0,"sec");

    // Rate : events/sec
    double rate = RATE;

    // Create observations
    GObservations obs;

    // Add some observation
    for (int i=0; i<6; ++i) {

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

        // Add events to the observation
        ob.events(events);
        ob.ontime(tmax.met()-tmin.met());
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
    GModelPar& result = ((obs.models())[0])[0];

    //check if converged
    test_assert(opt.status()==0, "Check if converged", "Optimizer did not convered"); 

    //check if value is correct
    test_value(result.value(),RATE,result.error()*3); 

    // Return
    return (((obs.models())[0])[0]);
}


/***********************************************************************//**
 * @brief Test optimizer with unbinned events and 1 thread
 ***************************************************************************/
void TestOpenMP::test_observations_optimizer_unbinned_1()
{
    // Test with 1 thread
    omp_set_num_threads(1);
    test_observations_optimizer(UN_BINNED);

    // Return
    return;
}

/***********************************************************************//**
 * @brief Test optimizer with unbinned events and 10 thread
 ***************************************************************************/
void TestOpenMP::test_observations_optimizer_unbinned_10()
{
    // Test with 10 threads
    omp_set_num_threads(10);
    test_observations_optimizer(UN_BINNED);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test optimizer with binned events and 1 thread
 ***************************************************************************/
void TestOpenMP::test_observations_optimizer_binned_1()
{
    // Test with 1 thread
    omp_set_num_threads(1);
    test_observations_optimizer(BINNED);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test optimizer with binned events and 10 threads
 ***************************************************************************/
void TestOpenMP::test_observations_optimizer_binned_10()
{
    // Test with 10 threads
    omp_set_num_threads(10);
    test_observations_optimizer(BINNED);

    // Return
    return;
}
#endif


/***************************************************************************
 * @brief Main entry point for test executable
 ***************************************************************************/
int main(void)
{
    // Allocate test suit container
    GTestSuites testsuites("Observation module");

    // Initially assume that we pass all tests
    bool success = true;

    // Create and append OpenMP test suite
    #ifdef _OPENMP
    TestOpenMP openmp;
    testsuites.append(openmp);
    #endif

    // Run the testsuites
    success = testsuites.run();

    // Save test report
    testsuites.save("reports/GObservation.xml");

    // Return success status
    return (success ? 0 : 1);
}
