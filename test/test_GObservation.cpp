/***************************************************************************
 *              test_GObservation.cpp - Test observation module            *
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
 * @author Jean-Baptiste Cayrou
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
//#include <iostream>
#include "testinst/GTestLib.hpp"
#include "test_GObservation.hpp"

/* __ OpenMP section _____________________________________________________ */
#ifdef _OPENMP
#include <omp.h>
#ifdef HAVE_OPENMP_DARWIN_KLUGE
#include <pthread.h>
pthread_attr_t gomp_thread_attr;
#endif
#endif

/* __ Coding definitions _________________________________________________ */
#define RATE      13.0        //!< Events per seconde. For events generation.
#define UN_BINNED 0
#define BINNED    1


/***********************************************************************//**
* @brief Set tests
***************************************************************************/
void TestGObservation::set(void)
{
    // Set test name
    name("Observation module");

    // Append tests
    append(static_cast<pfunction>(&TestGObservation::test_time_reference), "Test GTimeReference");
    append(static_cast<pfunction>(&TestGObservation::test_time), "Test GTime");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test GTimeReference
 ***************************************************************************/
void TestGObservation::test_time_reference(void)
{
    // Test void constructor
    test_try("Void constructor");
    try {
        GTimeReference reference;
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test copy constructor
    test_try("Copy constructor");
    try {
        GTimeReference reference;
        GTimeReference reference2(reference);
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test reference constructor
    test_try("Reference constructor");
    try {
        GTimeReference reference(55197.0, "s", "TT", "LOCAL");
        test_try_success();
        test_value(reference.mjdref(), 55197.0);
        test_assert(reference.timeunit() == "s",
                    "Time unit was \""+reference.timeunit()+"\", expected \"s\"");
        test_assert(reference.timesys() == "TT",
                    "Time system was \""+reference.timesys()+"\", expected \"TT\"");
        test_assert(reference.timeref() == "LOCAL",
                    "Time reference was \""+reference.timeref()+"\", expected \"LOCAL\"");
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test reference constructor
    test_try("Reference constructor (split reference)");
    try {
        GTimeReference reference(55197, 0.000766018518519, "s", "TT", "LOCAL");
        test_try_success();
        test_value(reference.mjdref(), 55197.000766018518519);
        test_assert(reference.timeunit() == "s",
                    "Time unit was \""+reference.timeunit()+"\", expected \"s\"");
        test_assert(reference.timesys() == "TT",
                    "Time system was \""+reference.timesys()+"\", expected \"TT\"");
        test_assert(reference.timeref() == "LOCAL",
                    "Time reference was \""+reference.timeref()+"\", expected \"LOCAL\"");
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test FITS file writing
    GTimeReference reference(55197.000766018518519, "s", "TT", "LOCAL");
    GFits          fits;
    GFitsBinTable  table;
    reference.write(&table);
    fits.append(table);
    fits.saveto("test_time_reference.fits", true);
    fits.close();

    // Read back from FITS file and check values
    fits.open("test_time_reference.fits");
    GFitsTable* hdu = fits.table(1);
    GTimeReference value(hdu);
    fits.close();
    test_value(value.mjdref(),  reference.mjdref());
    test_value(value.mjdrefi(), reference.mjdrefi());
    test_value(value.mjdreff(), reference.mjdreff());
    test_assert(value.timeunit() == reference.timeunit(),
                "Time unit was \""+value.timeunit()+"\", expected "+reference.timeunit()+".");
    test_assert(value.timesys() == reference.timesys(),
                "Time system was \""+value.timesys()+"\", expected "+reference.timesys()+".");
    test_assert(value.timeref() == reference.timeref(),
                "Time reference was \""+value.timeref()+"\", expected "+reference.timeref()+".");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test GTime
 ***************************************************************************/
void TestGObservation::test_time(void)
{
    // Test void constructor
    test_try("Void constructor");
    try {
        GTime time;
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test copy constructor
    test_try("Copy constructor");
    try {
        GTime time;
        GTime time2(time);
        test_try_success();
        test_assert(time == time2, "Time differs after using copy constructor.");
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test time constructor (seconds)
    test_try("Time constructor (seconds)");
    try {
        GTime time(1800.01);
        test_try_success();
        test_value(time.secs(), 1800.01);
        test_value(time.days(), 1800.01/86400.0);
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test time constructor (days)
    test_try("Time constructor (days)");
    try {
        GTime time(41.7, "days");
        test_try_success();
        test_value(time.days(), 41.7);
        test_value(time.secs(), 41.7*86400.0);
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test access methods
    double mjd_ref = 55197.000766018518519;
    double jd_ref  = 2455197.500766018518519;
    double t       = 123456.789;
    GTime time(t);
    test_value(time.jd(), t/86400.0 + jd_ref);
    test_value(time.mjd(), t/86400.0 + mjd_ref);
    test_value(time.secs(), t);
    test_value(time.days(), t/86400.0);

    // Test set method
    time.jd(57.9);
    test_value(time.jd(), 57.9);
    time.mjd(57.9);
    test_value(time.mjd(), 57.9);
    time.secs(57.9);
    test_value(time.secs(), 57.9);
    time.days(57.9);
    test_value(time.days(), 57.9);

    // Test convert method
    time.secs(t);
    test_value(time.convert(GTimeReference(55197.000766018518519, "days", "TT", "LOCAL")),
               t/86400.0);
    test_value(time.convert(GTimeReference(0.0, "s", "TT", "LOCAL")),
               t + mjd_ref*86400.0, 1.0e-6); //!< Poor precision on OpenSolaris

    // Test set method
    time.set(12.3, GTimeReference(55197.000766018518519, "days", "TT", "LOCAL"));
    test_value(time.days(), 12.3);
    time.set(12.3, GTimeReference(0.0, "secs", "TT", "LOCAL"));
    test_value(time.secs(), 12.3 - mjd_ref*86400.0);

    // Test operators
    GTime a(13.72);
    GTime b(6.28);
    test_value((a+b).secs(), 20.00);
    test_value((a-b).secs(), 7.44);
    test_value((a*3.3).secs(), 45.276);
    test_value((3.3*a).secs(), 45.276);
    test_value((a/13.72).secs(), 1.0);
    test_assert(a == a, "Equality operator corrupt.");
    test_assert(a != b, "Non-equality operator corrupt.");
    test_assert(a > b, "Greater than operator corrupt.");
    test_assert(a >= b, "Greater than or equal operator corrupt.");
    test_assert(b < a, "Less than operator corrupt.");
    test_assert(b <= a, "Less than or equal operator corrupt.");

    // Return
    return;
}


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
    GTime tmin(0.0);
    GTime tmax(1800.0);

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
    GModelPar& result = (*(obs.models()[0]))[0];

    //check if converged
    test_assert(opt.status()==0, "Check if converged", "Optimizer did not convered"); 

    //check if value is correct
    test_value(result.value(),RATE,result.error()*3); 

    // Return
    return (*(obs.models()[0]))[0];
}


/***********************************************************************//**
 * @brief Test optimizer with unbinned events and 1 thread
 ***************************************************************************/
void TestOpenMP::test_observations_optimizer_unbinned_1(void)
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
void TestOpenMP::test_observations_optimizer_unbinned_10(void)
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
void TestOpenMP::test_observations_optimizer_binned_1(void)
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
void TestOpenMP::test_observations_optimizer_binned_10(void)
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

    // Create and append test suites
    TestGObservation obs;
    testsuites.append(obs);

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
