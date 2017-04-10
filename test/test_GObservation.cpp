/***************************************************************************
 *              test_GObservation.cpp - Test observation module            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2017 by Jean-Baptiste Cayrou                        *
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
    append(static_cast<pfunction>(&TestGObservation::test_time_reference),
           "Test GTimeReference class");
    append(static_cast<pfunction>(&TestGObservation::test_time),
           "Test GTime class");
    append(static_cast<pfunction>(&TestGObservation::test_times),
           "Test GTimes class");
    append(static_cast<pfunction>(&TestGObservation::test_gti),
           "Test GGti class");
    append(static_cast<pfunction>(&TestGObservation::test_energy),
           "Test GEnergy class");
    append(static_cast<pfunction>(&TestGObservation::test_energies),
           "Test GEnergies class");
    append(static_cast<pfunction>(&TestGObservation::test_ebounds),
           "Test GEbounds class");
    append(static_cast<pfunction>(&TestGObservation::test_photons),
           "Test GPhotons class");
    append(static_cast<pfunction>(&TestGObservation::test_observations_optimizer),
           "Test GObservations::likelihood class");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone test suite
 *
 * @return Pointer to deep copy of test suite.
 ***************************************************************************/
TestGObservation* TestGObservation::clone(void) const
{
    // Clone test suite
    return new TestGObservation(*this);
}


/***********************************************************************//**
 * @brief Test GEbounds
 ***************************************************************************/
void TestGObservation::test_ebounds(void)
{
    // Test void constructor
    test_try("Void constructor");
    try {
        GEbounds ebds;
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Manipulate energy boudaries starting from an empty object
    GEbounds ebds;
    test_value(ebds.size(), 0, "GEbounds should have zero size.");
    test_assert(ebds.is_empty(), "GEbounds should be empty.");
    test_value(ebds.emin().MeV(), 0.0, 1.0e-10, "Minimum energy should be 0.");
    test_value(ebds.emax().MeV(), 0.0, 1.0e-10, "Maximum energy should be 0.");

    // Add empty interval
    ebds.append(GEnergy(1.0, "MeV"), GEnergy(1.0, "MeV"));
    test_value(ebds.size(), 1, "GEbounds should have 1 element.");
    test_assert(!ebds.is_empty(), "GEbounds should not be empty.");
    test_value(ebds.emin().MeV(), 1.0, 1.0e-10, "Minimum energy should be 1.");
    test_value(ebds.emax().MeV(), 1.0, 1.0e-10, "Maximum energy should be 1.");

    // Add one interval
    ebds.append(GEnergy(1.0, "MeV"), GEnergy(10.0, "MeV"));
    test_value(ebds.size(), 2, "GEbounds should have 2 elements.");
    test_assert(!ebds.is_empty(), "GEbounds should not be empty.");
    test_value(ebds.emin().MeV(), 1.0, 1.0e-10, "Minimum energy should be 1.");
    test_value(ebds.emax().MeV(), 10.0, 1.0e-10, "Maximum energy should be 10.");

    // Remove intervals
    ebds.remove(0);
    ebds.remove(0); // Now old intervals 1 is interval 0
    test_value(ebds.size(), 0, "GEbounds should have zero size.");
    test_assert(ebds.is_empty(), "GEbounds should be empty.");
    test_value(ebds.emin().MeV(), 0.0, 1.0e-10, "Minimum energy should be 0.");
    test_value(ebds.emax().MeV(), 0.0, 1.0e-10, "Maximum energy should be 0.");

    // Append two overlapping intervals
    ebds.append(GEnergy(1.0, "MeV"), GEnergy(100.0, "MeV"));
    ebds.append(GEnergy(10.0, "MeV"), GEnergy(1000.0, "MeV"));
    test_value(ebds.size(), 2, "GEbounds should have 2 elements.");
    test_assert(!ebds.is_empty(), "GEbounds should not be empty.");
    test_value(ebds.emin().MeV(), 1.0, 1.0e-10, "Minimum energy should be 1.");
    test_value(ebds.emax().MeV(), 1000.0, 1.0e-10, "Maximum energy should be 1000.");

    // Clear object
    ebds.clear();
    test_value(ebds.size(), 0, "GEbounds should have zero size.");
    test_assert(ebds.is_empty(), "GEbounds should be empty.");
    test_value(ebds.emin().MeV(), 0.0, 1.0e-10, "Minimum energy should be 0.");
    test_value(ebds.emax().MeV(), 0.0, 1.0e-10, "Maximum energy should be 0.");

    // Append two overlapping intervals in inverse order
    ebds.clear();
    ebds.append(GEnergy(10.0, "MeV"), GEnergy(1000.0, "MeV"));
    ebds.append(GEnergy(1.0, "MeV"), GEnergy(100.0, "MeV"));
    test_value(ebds.size(), 2, "GEbounds should have 2 elements.");
    test_assert(!ebds.is_empty(), "GEbounds should not be empty.");
    test_value(ebds.emin().MeV(), 1.0, 1.0e-10, "Minimum energy should be 1.");
    test_value(ebds.emax().MeV(), 1000.0, 1.0e-10, "Maximum energy should be 1000.");

    // Insert two overlapping intervals
    ebds.clear();
    ebds.insert(GEnergy(1.0, "MeV"), GEnergy(100.0, "MeV"));
    ebds.insert(GEnergy(10.0, "MeV"), GEnergy(1000.0, "MeV"));
    test_value(ebds.size(), 2, "GEbounds should have 2 elements.");
    test_assert(!ebds.is_empty(), "GEbounds should not be empty.");
    test_value(ebds.emin().MeV(), 1.0, 1.0e-10, "Minimum energy should be 1.");
    test_value(ebds.emax().MeV(), 1000.0, 1.0e-10, "Maximum energy should be 1000.");

    // Insert two overlapping intervals in inverse order
    ebds.clear();
    ebds.insert(GEnergy(10.0, "MeV"), GEnergy(1000.0, "MeV"));
    ebds.insert(GEnergy(1.0, "MeV"), GEnergy(100.0, "MeV"));
    test_value(ebds.size(), 2, "GEbounds should have 2 elements.");
    test_assert(!ebds.is_empty(), "GEbounds should not be empty.");
    test_value(ebds.emin().MeV(), 1.0, 1.0e-10, "Minimum energy should be 1.");
    test_value(ebds.emax().MeV(), 1000.0, 1.0e-10, "Maximum energy should be 1000.");

    // Merge two overlapping intervals
    ebds.clear();
    ebds.merge(GEnergy(1.0, "MeV"), GEnergy(100.0, "MeV"));
    ebds.merge(GEnergy(10.0, "MeV"), GEnergy(1000.0, "MeV"));
    test_value(ebds.size(), 1, "GEbounds should have 1 element.");
    test_assert(!ebds.is_empty(), "GEbounds should not be empty.");
    test_value(ebds.emin().MeV(), 1.0, 1.0e-10, "Minimum energy should be 1.");
    test_value(ebds.emax().MeV(), 1000.0, 1.0e-10, "Maximum energy should be 1000.");

    // Merge two overlapping intervals in inverse order
    ebds.clear();
    ebds.merge(GEnergy(10.0, "MeV"), GEnergy(1000.0, "MeV"));
    ebds.merge(GEnergy(1.0, "MeV"), GEnergy(100.0, "MeV"));
    test_value(ebds.size(), 1, "GEbounds should have 1 element.");
    test_assert(!ebds.is_empty(), "GEbounds should not be empty.");
    test_value(ebds.emin().MeV(), 1.0, 1.0e-10, "Minimum energy should be 1.");
    test_value(ebds.emax().MeV(), 1000.0, 1.0e-10, "Maximum energy should be 1000.");

    // Check linear boundaries
    ebds.clear();
    ebds.set_lin(2, GEnergy(1.0, "MeV"), GEnergy(3.0, "MeV"));
    test_value(ebds.size(), 2, "GEbounds should have 2 elements.");
    test_assert(!ebds.is_empty(), "GEbounds should not be empty.");
    test_value(ebds.emin(0).MeV(), 1.0, 1.0e-10, "Bin 0 minimum energy should be 1.");
    test_value(ebds.emin(1).MeV(), 2.0, 1.0e-10, "Bin 1 minimum energy should be 2.");
    test_value(ebds.emax(0).MeV(), 2.0, 1.0e-10, "Bin 0 maximum energy should be 2.");
    test_value(ebds.emax(1).MeV(), 3.0, 1.0e-10, "Bin 1 maximum energy should be 3.");
    test_value(ebds.emin().MeV(), 1.0, 1.0e-10, "Minimum energy should be 1.");
    test_value(ebds.emax().MeV(), 3.0, 1.0e-10, "Maximum energy should be 3.");

    // Check logarithmic boundaries
    ebds.clear();
    ebds.set_log(2, GEnergy(1.0, "MeV"), GEnergy(100.0, "MeV"));
    test_value(ebds.size(), 2, "GEbounds should have 2 elements.");
    test_assert(!ebds.is_empty(), "GEbounds should not be empty.");
    test_value(ebds.emin(0).MeV(), 1.0, 1.0e-10, "Bin 0 minimum energy should be 1.");
    test_value(ebds.emin(1).MeV(), 10.0, 1.0e-10, "Bin 1 minimum energy should be 10.");
    test_value(ebds.emax(0).MeV(), 10.0, 1.0e-10, "Bin 0 maximum energy should be 10.");
    test_value(ebds.emax(1).MeV(), 100.0, 1.0e-10, "Bin 1 maximum energy should be 100.");
    test_value(ebds.emin().MeV(), 1.0, 1.0e-10, "Minimum energy should be 1.");
    test_value(ebds.emax().MeV(), 100.0, 1.0e-10, "Maximum energy should be 100.");

    // Check boundary extension
    GEbounds ext(1, GEnergy(100.0, "MeV"), GEnergy(1000.0, "MeV"));
    ebds.extend(ext);
    test_value(ebds.size(), 3, "GEbounds should have 3 elements.");
    test_assert(!ebds.is_empty(), "GEbounds should not be empty.");
    test_value(ebds.emin(0).MeV(), 1.0, 1.0e-10, "Bin 0 minimum energy should be 1.");
    test_value(ebds.emin(1).MeV(), 10.0, 1.0e-10, "Bin 1 minimum energy should be 10.");
    test_value(ebds.emin(2).MeV(), 100.0, 1.0e-10, "Bin 2 minimum energy should be 100.");
    test_value(ebds.emax(0).MeV(), 10.0, 1.0e-10, "Bin 0 maximum energy should be 10.");
    test_value(ebds.emax(1).MeV(), 100.0, 1.0e-10, "Bin 1 maximum energy should be 100.");
    test_value(ebds.emax(2).MeV(), 1000.0, 1.0e-10, "Bin 1 maximum energy should be 1000.");
    test_value(ebds.emin().MeV(), 1.0, 1.0e-10, "Minimum energy should be 1.");
    test_value(ebds.emax().MeV(), 1000.0, 1.0e-10, "Maximum energy should be 1000.");

    // Check emean, elogmean and ewidth methods
    ebds.set_log(1, GEnergy(1.0, "MeV"), GEnergy(10.0, "MeV"));
    test_value(ebds.emean(0).MeV(), 5.5, 1.0e-10, "Mean energy should be 5.5.");
    test_value(ebds.elogmean(0).MeV(), 3.16227766017, 1.0e-10, "Log mean energy should be 3.16227766017.");
    test_value(ebds.ewidth(0).MeV(), 9.0, 1.0e-10, "Energy width should be 9.0.");

    // Check appending of invalid interval
    test_try("Test appending of invalid interval");
    try {
        ebds.append(GEnergy(100.0, "MeV"), GEnergy(10.0, "MeV"));
        test_try_failure("Appending an invalid interval shall throw an exception.");
    }
    catch (GException::invalid_argument &e) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Check constructing of linear invalid energy boundaries
    test_try("Test constructing of linear invalid energy boundaries");
    try {
        GEbounds bad(10, GEnergy(100.0, "MeV"), GEnergy(10.0, "MeV"), false);
        test_try_failure("Constructing of linear invalid energy boundaries "
                         "shall throw an exception.");
    }
    catch (GException::invalid_argument &e) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Check constructing of logarithmic invalid energy boundaries
    test_try("Test constructing of logarithmic invalid energy boundaries");
    try {
        GEbounds bad(10, GEnergy(100.0, "MeV"), GEnergy(10.0, "MeV"), true);
        test_try_failure("Constructing of logarithmic invalid energy "
                         "boundaries shall throw an exception.");
    }
    catch (GException::invalid_argument &e) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Check constructing of logarithmic invalid energy boundaries
    test_try("Test constructing of logarithmic invalid energy boundaries");
    try {
        GEbounds bad(10, GEnergy(0.0, "MeV"), GEnergy(10.0, "MeV"), true);
        test_try_failure("Constructing of logarithmic invalid energy "
                         "boundaries shall throw an exception.");
    }
    catch (GException::invalid_argument &e) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Check constructing of logarithmic invalid energy boundaries
    test_try("Test constructing of logarithmic invalid energy boundaries");
    try {
        GEbounds bad(10, GEnergy(10.0, "MeV"), GEnergy(0.0, "MeV"), true);
        test_try_failure("Constructing of logarithmic invalid energy "
                         "boundaries shall throw an exception.");
    }
    catch (GException::invalid_argument &e) {
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Check energy containment
    GEbounds containment(1, GEnergy(100.0, "MeV"), GEnergy(1000.0, "MeV"));
    test_assert(containment.contains(GEnergy(200.0, "MeV")),
                "Energy 200 MeV should be contained.");
    test_assert(!containment.contains(GEnergy(10.0, "MeV")),
                "Energy 10 MeV should not be contained.");
    test_assert(!containment.contains(GEnergy(2000.0, "MeV")),
                "Energy 2000 MeV should not be contained.");
    test_assert(containment.contains(GEnergy(200.0, "MeV"), GEnergy(800.0, "MeV")),
                "Energy bin [200,800] MeV should be contained.");
    test_assert(!containment.contains(GEnergy(80.0, "MeV"), GEnergy(800.0, "MeV")),
                "Energy bin [80,800] MeV should not be contained.");
    test_assert(!containment.contains(GEnergy(200.0, "MeV"), GEnergy(2000.0, "MeV")),
                "Energy bin [200,2000] MeV should not be contained.");
    test_assert(!containment.contains(GEnergy(80.0, "MeV"), GEnergy(2000.0, "MeV")),
                "Energy bin [80,2000] MeV should not be contained.");

    // Remove test file
    GFilename filename("test_ebounds.fits");
    filename.remove();

    // Check saving
    containment.save("test_ebounds.fits");
    GEbounds load1("test_ebounds.fits");
    test_value(load1.size(), 1, "GEbounds should have 1 element.");
    test_value(load1.emin().MeV(), 100.0, 1.0e-10, "Minimum energy should be 100.");
    test_value(load1.emax().MeV(), 1000.0, 1.0e-10, "Maximum energy should be 1000.");

    // Check saving in a different extnsion
    ebds.save("test_ebounds.fits[ENERGY BOUNDARIES 2]", true);
    GEbounds load2("test_ebounds.fits[ENERGY BOUNDARIES 2]");
    test_value(load2.size(), 1, "GEbounds should have 1 element.");
    test_value(load2.emin().MeV(), 1.0, 1.0e-10, "Minimum energy should be 1.");
    test_value(load2.emax().MeV(), 10.0, 1.0e-10, "Maximum energy should be 10.");

    // Check XML write and read methods
    GXmlElement element;
    ebds.write(element);
    GEbounds xml1(element);
    test_value(xml1.size(), 1, "XML write and read methods");
    test_value(xml1.emin().MeV(),  1.0, 1.0e-10, "XML write and read methods");
    test_value(xml1.emax().MeV(), 10.0, 1.0e-10, "XML write and read methods");

    // Check energies constructor method
    GEnergies energies1(3, GEnergy(1.0, "MeV"), GEnergy(10.0, "MeV"), false);
    GEbounds ebds1(energies1);
    test_value(ebds1.size(), 2, "GEnergies constructor (3 elements)");
    test_value(ebds1.emin().MeV(),  1.0, 1.0e-10, "GEnergies constructor (3 elements)");
    test_value(ebds1.emax(0).MeV(), 5.5, 1.0e-10, "GEnergies constructor (3 elements)");
    test_value(ebds1.emax().MeV(), 10.0, 1.0e-10, "GEnergies constructor (3 elements)");

    // Check energies constructor method
    GEnergies energies2(1, GEnergy(1.0, "MeV"), GEnergy(1.0, "MeV"), false);
    GEbounds ebds2(energies2);
    test_value(ebds2.size(), 1, "GEnergies constructor (1 element)");
    test_value(ebds2.emin().MeV(), 1.0, 1.0e-10, "GEnergies constructor (1 element)");
    test_value(ebds2.emax().MeV(), 1.0, 1.0e-10, "GEnergies constructor (1 element)");

    // Check energy boundary setter method
    ebds2.emin(0, GEnergy(123.0, "MeV"));
    ebds2.emax(0, GEnergy(511.0, "MeV"));
    test_value(ebds2.size(), 1, "Energy boundary setter");
    test_value(ebds2.emin().MeV(), 123.0, 1.0e-10, "Energy boundary setter");
    test_value(ebds2.emax().MeV(), 511.0, 1.0e-10, "Energy boundary setter");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test GGti
 ***************************************************************************/
void TestGObservation::test_gti(void)
{
    // Test void constructor
    test_try("Void constructor");
    try {
        GGti gti;
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Manipulate GTIs starting from an empty object
    GGti gti;
    test_value(gti.size(), 0, "GGti should have zero size.");
    test_assert(gti.is_empty(), "GGti should be empty.");
    test_value(gti.tstart().secs(), 0.0, 1.0e-7, "Start time should be 0.");
    test_value(gti.tstop().secs(), 0.0, 1.0e-7, "Stop time should be 0.");

    // Add empty interval
    gti.append(GTime(1.0), GTime(1.0));
    test_value(gti.size(), 1, "GGti should have 1 interval.");
    test_assert(!gti.is_empty(), "GGti should not be empty.");
    test_value(gti.tstart().secs(), 1.0, 1.0e-7, "Start time should be 1.");
    test_value(gti.tstop().secs(), 1.0, 1.0e-7, "Stop time should be 1.");

    // Add one interval
    gti.append(GTime(1.0), GTime(10.0));
    test_value(gti.size(), 2, "GGti should have 2 intervals.");
    test_assert(!gti.is_empty(), "GGti should not be empty.");
    test_value(gti.tstart().secs(), 1.0, 1.0e-7, "Start time should be 1.");
    test_value(gti.tstop().secs(), 10.0, 1.0e-7, "Stop time should be 10.");

    // Remove interval
    gti.remove(0);
    gti.remove(0);
    test_value(gti.size(), 0, "GGti should have zero size.");
    test_assert(gti.is_empty(), "GGti should be empty.");
    test_value(gti.tstart().secs(), 0.0, 1.0e-7, "Start time should be 0.");
    test_value(gti.tstop().secs(), 0.0, 1.0e-7, "Stop time should be 0.");

    // Append two overlapping intervals
    gti.append(GTime(1.0), GTime(100.0));
    gti.append(GTime(10.0), GTime(1000.0));
    test_value(gti.size(), 2, "GGti should have 2 intervals.");
    test_assert(!gti.is_empty(), "GGti should not be empty.");
    test_value(gti.tstart().secs(), 1.0, 1.0e-7, "Start time should be 1.");
    test_value(gti.tstop().secs(), 1000.0, 1.0e-7, "Stop time should be 1000.");

    // Clear object
    gti.clear();
    test_value(gti.size(), 0, "GGti should have zero size.");
    test_assert(gti.is_empty(), "GGti should be empty.");
    test_value(gti.tstart().secs(), 0.0, 1.0e-7, "Start time should be 0.");
    test_value(gti.tstop().secs(), 0.0, 1.0e-7, "Stop time should be 0.");

    // Append two overlapping intervals in inverse order
    gti.clear();
    gti.append(GTime(10.0), GTime(1000.0));
    gti.append(GTime(1.0), GTime(100.0));
    test_value(gti.size(), 2, "GGti should have 2 intervals.");
    test_assert(!gti.is_empty(), "GGti should not be empty.");
    test_value(gti.tstart().secs(), 1.0, 1.0e-7, "Start time should be 1.");
    test_value(gti.tstop().secs(), 1000.0, 1.0e-7, "Stop time should be 1000.");

    // Insert two overlapping intervals
    gti.clear();
    gti.insert(GTime(1.0), GTime(100.0));
    gti.insert(GTime(10.0), GTime(1000.0));
    test_value(gti.size(), 2, "GGti should have 2 intervals.");
    test_assert(!gti.is_empty(), "GGti should not be empty.");
    test_value(gti.tstart().secs(), 1.0, 1.0e-10, "Start time should be 1.");
    test_value(gti.tstop().secs(), 1000.0, 1.0e-10, "Stop time should be 1000.");

    // Insert two overlapping intervals in inverse order
    gti.clear();
    gti.insert(GTime(10.0), GTime(1000.0));
    gti.insert(GTime(1.0), GTime(100.0));
    test_value(gti.size(), 2, "GGti should have 2 intervals.");
    test_assert(!gti.is_empty(), "GGti should not be empty.");
    test_value(gti.tstart().secs(), 1.0, 1.0e-7, "Start time should be 1.");
    test_value(gti.tstop().secs(), 1000.0, 1.0e-7, "Stop time should be 1000.");

    // Merge two overlapping intervals
    gti.clear();
    gti.merge(GTime(1.0), GTime(100.0));
    gti.merge(GTime(10.0), GTime(1000.0));
    test_value(gti.size(), 1, "GGti should have 1 interval.");
    test_assert(!gti.is_empty(), "GGti should not be empty.");
    test_value(gti.tstart().secs(), 1.0, 1.0e-7, "Start time should be 1.");
    test_value(gti.tstop().secs(), 1000.0, 1.0e-7, "Stop time should be 1000.");

    // Merge two overlapping intervals in inverse order
    gti.clear();
    gti.merge(GTime(10.0), GTime(1000.0));
    gti.merge(GTime(1.0), GTime(100.0));
    test_value(gti.size(), 1, "GGti should have 1 interval.");
    test_assert(!gti.is_empty(), "GGti should not be empty.");
    test_value(gti.tstart().secs(), 1.0, 1.0e-7, "Start time should be 1.");
    test_value(gti.tstop().secs(), 1000.0, 1.0e-7, "Stop time should be 1000.");

    // Check extension
    gti.clear();
    gti.append(GTime(1.0), GTime(10.0));
    gti.append(GTime(10.0), GTime(100.0));
    GGti ext;
    ext.append(GTime(100.0), GTime(1000.0));
    gti.extend(ext);
    test_value(gti.size(), 3, "GGti should have 3 intervals.");
    test_assert(!gti.is_empty(), "GGti should not be empty.");
    test_value(gti.tstart(0).secs(), 1.0, 1.0e-7, "Bin 0 start time should be 1.");
    test_value(gti.tstart(1).secs(), 10.0, 1.0e-7, "Bin 1 start time should be 10.");
    test_value(gti.tstart(2).secs(), 100.0, 1.0e-7, "Bin 2 start time should be 100.");
    test_value(gti.tstop(0).secs(), 10.0, 1.0e-7, "Bin 0 stop time should be 10.");
    test_value(gti.tstop(1).secs(), 100.0, 1.0e-7, "Bin 1 stop time should be 100.");
    test_value(gti.tstop(2).secs(), 1000.0, 1.0e-7, "Bin 2 stop time should be 1000.");
    test_value(gti.tstart().secs(), 1.0, 1.0e-7, "Start time should be 1.");
    test_value(gti.tstop().secs(), 1000.0, 1.0e-7, "Stop time should be 1000.");

    // Remove test file
    GFilename file1("test_gti1.fits");
    GFilename file2("test_gti2.fits");
    file1.remove();
    file2.remove();

    // Check saving in and loading from FITS file
    gti.save("test_gti1.fits", true);
    GGti test1("test_gti1.fits");
    test_value(test1.size(), 3, "GGti should have 3 intervals.");
    test_value(test1.tstart(0).secs(), 1.0, 1.0e-7, "Bin 0 start time should be 1.");
    test_value(test1.tstart(1).secs(), 10.0, 1.0e-7, "Bin 1 start time should be 10.");
    test_value(test1.tstart(2).secs(), 100.0, 1.0e-7, "Bin 2 start time should be 100.");
    test_value(test1.tstop(0).secs(), 10.0, 1.0e-7, "Bin 0 stop time should be 10.");
    test_value(test1.tstop(1).secs(), 100.0, 1.0e-7, "Bin 1 stop time should be 100.");
    test_value(test1.tstop(2).secs(), 1000.0, 1.0e-7, "Bin 2 stop time should be 1000.");
    test_value(test1.tstart().secs(), 1.0, 1.0e-7, "Start time should be 1.");
    test_value(test1.tstop().secs(), 1000.0, 1.0e-7, "Stop time should be 1000.");

    // Check saving in and loading from FITS file with specific extension
    gti.save("test_gti2.fits[GTI_2]", true);
    GGti test2("test_gti2.fits[GTI_2]");
    test_value(test2.size(), 3, "GGti should have 3 intervals.");
    test_value(test2.tstart(0).secs(), 1.0, 1.0e-7, "Bin 0 start time should be 1.");
    test_value(test2.tstart(1).secs(), 10.0, 1.0e-7, "Bin 1 start time should be 10.");
    test_value(test2.tstart(2).secs(), 100.0, 1.0e-7, "Bin 2 start time should be 100.");
    test_value(test2.tstop(0).secs(), 10.0, 1.0e-7, "Bin 0 stop time should be 10.");
    test_value(test2.tstop(1).secs(), 100.0, 1.0e-7, "Bin 1 stop time should be 100.");
    test_value(test2.tstop(2).secs(), 1000.0, 1.0e-7, "Bin 2 stop time should be 1000.");
    test_value(test2.tstart().secs(), 1.0, 1.0e-7, "Start time should be 1.");
    test_value(test2.tstop().secs(), 1000.0, 1.0e-7, "Stop time should be 1000.");

    // Check reading from XML element using "tmin" and "tmax" format
    GXmlElement element("observation");
    element.append(GXmlElement("parameter name=\"GoodTimeIntervals\" tmin=\"0\" tmax=\"100.0\""));
    element.append(GXmlElement("parameter name=\"TimeReference\" mjdrefi=\"51544\" mjdreff=\"0.5\" timeunit=\"s\" timesys=\"TT\" timeref=\"LOCAL\""));
    GGti test3(element);
    test_value(test3.size(), 1, "GGti should have 1 interval.");
    test_value(test3.tstart(0).convert(GTimeReference(51544.5, "s")), 0.0, 1.0e-7, "Bin 0 start time should be 0.");
    test_value(test3.tstop(0).convert(GTimeReference(51544.5, "s")), 100.0, 1.0e-7, "Bin 0 stop time should be 100.");
    test_value(test3.tstart().convert(GTimeReference(51544.5, "s")), 0.0, 1.0e-7, "Start time should be 0.");
    test_value(test3.tstop().convert(GTimeReference(51544.5, "s")), 100.0, 1.0e-7, "Stop time should be 100.");

    // Check writing of GTI into XML element
    element.clear();
    test3.write(element);
    const GXmlElement* par = gammalib::xml_get_par("test_GObservation", element, "GoodTimeIntervals");
    test_value(gammalib::tofloat(par->attribute("tmin")), 0.0, 1.0e-7,
               "Attribute tmin should be 0, found "+
               par->attribute("tmin")+".");
    test_value(gammalib::tofloat(par->attribute("tmax")), 100.0, 1.0e-7,
               "Attribute tmax should be 100, found "+
               par->attribute("tmax")+".");
    par = gammalib::xml_get_par("test_GObservation", element, "TimeReference");
    test_value(gammalib::tofloat(par->attribute("mjdrefi")), 51544.0, 1.0e-7,
               "Attribute mjdrefi should be 51544, found "+
               par->attribute("mjdrefi")+".");
    test_value(gammalib::tofloat(par->attribute("mjdreff")), 0.5, 1.0e-7,
               "Attribute mjdreff should be 0.5, found "+
               par->attribute("mjdreff")+".");
    test_assert(par->attribute("timeunit") == "s",
               "Attribute timeunit should be \"s\", found "+
               par->attribute("timeunit")+".");
    test_assert(par->attribute("timesys") == "TT",
               "Attribute timesys should be \"TT\", found "+
               par->attribute("timesys")+".");
    test_assert(par->attribute("timeref") == "LOCAL",
               "Attribute timeref should be \"LOCAL\", found "+
               par->attribute("timeref")+".");

    // Check reading from XML element using "file" format
    element = GXmlElement("observation");
    element.append(GXmlElement("parameter name=\"GoodTimeIntervals\" file=\"test_gti2.fits[GTI_2]\""));
    GGti test4(element);
    test_value(test4.size(), 3, "GGti should have 3 intervals.");
    test_value(test4.tstart(0).secs(), 1.0, 1.0e-7, "Bin 0 start time should be 1.");
    test_value(test4.tstart(1).secs(), 10.0, 1.0e-7, "Bin 1 start time should be 10.");
    test_value(test4.tstart(2).secs(), 100.0, 1.0e-7, "Bin 2 start time should be 100.");
    test_value(test4.tstop(0).secs(), 10.0, 1.0e-7, "Bin 0 stop time should be 10.");
    test_value(test4.tstop(1).secs(), 100.0, 1.0e-7, "Bin 1 stop time should be 100.");
    test_value(test4.tstop(2).secs(), 1000.0, 1.0e-7, "Bin 2 stop time should be 1000.");
    test_value(test4.tstart().secs(), 1.0, 1.0e-7, "Start time should be 1.");
    test_value(test4.tstop().secs(), 1000.0, 1.0e-7, "Stop time should be 1000.");

    // Append one time interval to test change of GTI file
    test4.append(GTime(1000.0), GTime(10000.0));

    // Check writing of GTI into XML element
    element.clear();
    test4.write(element);
    par = gammalib::xml_get_par("test_GObservation", element, "GoodTimeIntervals");
    test_assert(par->attribute("file") == "test_gti2.fits[GTI_2]",
               "Attribute file should be \"test_gti2.fits[GTI_2]\", found "+
               par->attribute("file")+".");
    GGti test5(par->attribute("file"));
    test_value(test5.size(), 4, "GGti should have 4 intervals.");
    test_value(test5.tstart(0).secs(), 1.0, 1.0e-7, "Bin 0 start time should be 1.");
    test_value(test5.tstart(1).secs(), 10.0, 1.0e-7, "Bin 1 start time should be 10.");
    test_value(test5.tstart(2).secs(), 100.0, 1.0e-7, "Bin 2 start time should be 100.");
    test_value(test5.tstart(3).secs(), 1000.0, 1.0e-7, "Bin 3 start time should be 1000.");
    test_value(test5.tstop(0).secs(), 10.0, 1.0e-7, "Bin 0 stop time should be 10.");
    test_value(test5.tstop(1).secs(), 100.0, 1.0e-7, "Bin 1 stop time should be 100.");
    test_value(test5.tstop(2).secs(), 1000.0, 1.0e-7, "Bin 2 stop time should be 1000.");
    test_value(test5.tstop(3).secs(), 10000.0, 1.0e-7, "Bin 3 stop time should be 10000.");
    test_value(test5.tstart().secs(), 1.0, 1.0e-7, "Start time should be 1.");
    test_value(test5.tstop().secs(), 10000.0, 1.0e-7, "Stop time should be 10000.");

    // Remove test file
    GFilename filename("test_gti.fits");
    filename.remove();

    // Check saving
    test4.save("test_gti.fits");
    GGti load1("test_gti.fits");
    test_value(load1.size(), 4, "GGti should have 4 intervals.");
    test_value(load1.tstart().secs(), 1.0, 1.0e-7, "Start time should be 1.");
    test_value(load1.tstop().secs(), 10000.0, 1.0e-7, "Stop time should be 10000.");

    // Check saving in a different extnsion
    test3.save("test_gti.fits[GOOD TIME INTERVALS]", true);
    GGti load2("test_gti.fits[GOOD TIME INTERVALS]");
    test_value(load2.size(), 1, "GGti should have 1 interval.");
    test_value(load2.tstart().convert(GTimeReference(51544.5, "s")), 0.0, 1.0e-7, "Start time should be 0.");
    test_value(load2.tstop().convert(GTimeReference(51544.5, "s")), 100.0, 1.0e-7, "Stop time should be 100.");

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
    reference.write(table);
    fits.append(table);
    fits.saveto("test_time_reference.fits", true);
    fits.close();

    // Read back from FITS file and check values
    fits.open("test_time_reference.fits");
    const GFitsTable& hdu = *fits.table(1);
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

    // Test time access methods
    // Reference times are from
    // - http://heasarc.gsfc.nasa.gov/cgi-bin/Tools/xTime/xTime.pl
    GTime time(123456.789);
    test_value(time.jd(), 2455198.92966648);
    test_value(time.jd("TT"), 2455198.92966648);
    test_value(time.jd("TAI"), 2455198.92929398);   // TAI=TT-0.0003725
    test_value(time.jd("UTC"), 2455198.92890046);
    test_value(time.mjd(), 55198.42966648);
    test_value(time.mjd("TT"), 55198.42966648);
    test_value(time.mjd("TAI"), 55198.42929398);    // TAI=TT-0.0003725
    test_value(time.mjd("UTC"), 55198.42890046);
    test_value(time.secs(), 123456.789);
    test_value(time.secs("TT"), 123456.789);
    test_value(time.secs("TAI"), 123424.605);       // TAI=TT-32.184
    test_value(time.secs("UTC"), 123390.60487);     // UTC=TT-66.184126
    test_value(time.days(), 1.4288979);
    test_value(time.days("TT"), 1.4288979);
    test_value(time.days("TAI"), 1.4285254);        // TAI=TT-0.0003725
    test_value(time.days("UTC"), 1.42813188);       // UTC=TT-0.00076602
    test_assert(time.utc() == "2010-01-02T10:17:37",
                "GTime::utc(): 2010-01-02T10:17:37 expected, "+
                time.utc()+" found.");

    // Test set method. The UTC conversion to MJD (TT) was done using xTime.
    // Note that xTime missed a second for the "2006-01-01T00:00:00"
    // conversion, hence I added the second by hand to the expected
    // result
    time.jd(57.9);
    test_value(time.jd(), 57.9);
    time.mjd(57.9);
    test_value(time.mjd(), 57.9);
    time.secs(57.9);
    test_value(time.secs(), 57.9);
    time.days(57.9);
    test_value(time.days(), 57.9);
    time.utc("1994-01-01T00:00:00");
    test_value(time.utc(), "1994-01-01T00:00:00");
    test_value(time.mjd(), 49353.00069657, 1.0e-6);
    time.utc("1999-12-31T23:59:59");
    test_value(time.utc(), "1999-12-31T23:59:59");
    test_value(time.mjd(), 51544.00073130, 1.0e-6);
    time.utc("2000-01-01T00:00:00");
    test_value(time.utc(), "2000-01-01T00:00:00");
    test_value(time.mjd(), 51544.00074287, 1.0e-6);
    time.utc("2006-01-01T00:00:00");
    test_value(time.utc(), "2006-01-01T00:00:00");
    test_value(time.mjd(), 53736.0007544, 1.0e-6);
    time.utc("2014-10-12T22:08:37");
    test_value(time.utc(), "2014-10-12T22:08:37");
    test_value(time.mjd(), 56942.92342806, 1.0e-6);

    // Test Julian Day set method for various time systems
    time.jd(2455198.92966648, "TT");
    test_value(time.jd(), 2455198.92966648);
    time.jd(2455198.92929398, "TAI");
    test_value(time.jd(), 2455198.92966648);
    time.jd(2455198.92890046, "UTC");
    test_value(time.jd(), 2455198.92966648);

    // Test Modified Julian Day set method for various time systems
    time.mjd(55198.42966648, "TT");
    test_value(time.mjd(), 55198.42966648);
    time.mjd(55198.42929398, "TAI");
    test_value(time.mjd(), 55198.42966648);
    time.mjd(55198.42890046, "UTC");
    test_value(time.mjd(), 55198.42966648);

    // Test seconds set method for various time systems
    time.secs(123456.789, "TT");
    test_value(time.secs(), 123456.789);
    time.secs(123424.605, "TAI");
    test_value(time.secs(), 123456.789);
    time.secs(123390.60487, "UTC");
    test_value(time.secs(), 123456.789);

    // Test days set method for various time systems
    time.days(1.4288979, "TT");
    test_value(time.days(), 1.4288979);
    time.days(1.4285254, "TAI");
    test_value(time.days(), 1.4288979);
    time.days(1.42813188, "UTC");
    test_value(time.days(), 1.4288979);

    // Test convert method
    double mjd_ref = 55197.000766018518519;
    time.secs(123456.789);
    test_value(time.convert(GTimeReference(mjd_ref, "days", "TT", "LOCAL")),
               1.4288979);
    test_value(time.convert(GTimeReference(mjd_ref, "days", "TAI", "LOCAL")),
               1.4285254);
    test_value(time.convert(GTimeReference(mjd_ref, "days", "UTC", "LOCAL")),
               1.42813188);
    test_value(time.convert(GTimeReference(0.0, "s", "TT", "LOCAL")),
               123456.789 + mjd_ref*86400.0, 1.0e-6); //!< Poor precision on OpenSolaris

    // Test set method
    time.set(1.4288979, GTimeReference(mjd_ref, "days", "TT", "LOCAL"));
    test_value(time.days(), 1.4288979);
    time.set(1.4285254, GTimeReference(mjd_ref, "days", "TAI", "LOCAL"));
    test_value(time.days(), 1.4288979);
    time.set(1.42813188, GTimeReference(mjd_ref, "days", "UTC", "LOCAL"));
    test_value(time.days(), 1.4288979);
    time.set(12.3, GTimeReference(0.0, "secs", "TT", "LOCAL"));
    test_value(time.secs(), 12.3 - mjd_ref*86400.0);

    // Test conversion to different time systems
    time.utc("2005-10-08T14:30:25");
    test_value(time.convert(GTimeReference(0.0, "days", "TT", "LOCAL")),
               53651.60519889, 1.0e-6);
    test_value(time.convert(GTimeReference(0.0, "days", "UTC", "LOCAL")),
               53651.60445602, 1.0e-6);

    // Test setting to different time systems
    time.set(53651.60519889, GTimeReference(0.0, "days", "TT", "LOCAL"));
    test_value(time.utc(), "2005-10-08T14:30:25");
    time.set(53651.60445602, GTimeReference(0.0, "days", "UTC", "LOCAL"));
    test_value(time.utc(), "2005-10-08T14:30:25");

    // Test now method (just test that it does not core dump; cannot really
    // test the value :-)
    time.now();

    // Test gmst and gast methods. Precision is 1 second since the methods
    // use UTC instead of UT1. Therefore  use a precision of 2.777e-4 in
    // the unit tests. The actual precision is even better (not sure why).
    // References from http://dc.zah.uni-heidelberg.de/apfs/times/q/form
    time.utc("1994-07-02T09:10:11");     // Largest UTC-UT1 difference
    test_value(time.gmst(), 3.8483111, 2.777e-4);
    test_value(time.gast(), 3.8485481, 2.777e-4);
    test_value(time.lmst(123.456), 19.617911, 2.777e-4);
    test_value(time.last(123.456), 19.618148, 2.777e-4);
    time.utc("1997-11-11T11:11:11");
    test_value(time.gmst(), 14.562163, 2.777e-4);
    test_value(time.gast(), 14.562068, 2.777e-4);
    test_value(time.lmst(123.456), 6.331763, 2.777e-4);
    test_value(time.last(123.456), 6.331668, 2.777e-4);
    time.utc("2008-10-04T10:30:23");
    test_value(time.gmst(), 11.405402, 2.777e-4);
    test_value(time.gast(), 11.405589, 2.777e-4);
    test_value(time.lmst(123.456), 3.175002, 2.777e-4);
    test_value(time.last(123.456), 3.175189, 2.777e-4);
    time.utc("2016-02-29T23:59:59");
    test_value(time.gmst(), 10.615045, 2.777e-4);
    test_value(time.gast(), 10.615018, 2.777e-4);
    test_value(time.lmst(123.456), 2.384645, 2.777e-4);
    test_value(time.last(123.456), 2.384618, 2.777e-4);

    // Test operators
    GTime a(13.72);
    GTime b(6.28);
    test_value((a+6.28).secs(), 20.00, 1.0e-6, "Seconds right addition operator");
    test_value((6.28+a).secs(), 20.00, 1.0e-6, "Seconds left addition operator");
    test_value((a-6.28).secs(), 7.44, 1.0e-6, "Seconds substraction operator");
    test_value((a-b), 7.44, 1.0e-6, "Time subtraction operator");
    GTime c = a;
    c += 6.28,
    test_value(c.secs(), 20.00, 1.0e-6, "Seconds unary addition operator");
    GTime d = a;
    d -= 6.28,
    test_value(d.secs(), 7.44, 1.0e-6, "Seconds unary subtraction operator");
    test_assert(a == a, "Equality operator");
    test_assert(a != b, "Non-equality operator");
    test_assert(a > b, "Greater than operator");
    test_assert(a >= b, "Greater than or equal operator");
    test_assert(b < a, "Less than operator");
    test_assert(b <= a, "Less than or equal operator");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test GTimes
 ***************************************************************************/
void TestGObservation::test_times(void)
{
    // Test void constructor
    test_try("Void constructor");
    try {
        GTimes times;
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Manipulate GTimes starting from an empty object
    GTimes times;
    test_value(times.size(), 0, "GTimes should have zero size.");
    test_assert(times.is_empty(), "GTimes should be empty.");

    // Add a time
    times.append(GTime());
    test_value(times.size(), 1, "GTimes should have 1 time.");
    test_assert(!times.is_empty(), "GTimes should not be empty.");

    // Remove time
    times.remove(0);
    test_value(times.size(), 0, "GTimes should have zero size.");
    test_assert(times.is_empty(), "GTimes should be empty.");

    // Append two times
    times.append(GTime());
    times.append(GTime());
    test_value(times.size(), 2, "GTimes should have 2 times.");
    test_assert(!times.is_empty(), "GTimes should not be empty.");

    // Clear object
    times.clear();
    test_value(times.size(), 0, "GTimes should have zero size.");
    test_assert(times.is_empty(), "GTimes should be empty.");

    // Insert two times
    times.insert(0, GTime());
    times.insert(0, GTime());
    test_value(times.size(), 2, "GTimes should have 2 times.");
    test_assert(!times.is_empty(), "GTimes should not be empty.");

    // Extend times
    times.extend(times);
    test_value(times.size(), 4, "GTimes should have 4 times.");
    test_assert(!times.is_empty(), "GTimes should not be empty.");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test GEnergy class
 ***************************************************************************/
void TestGObservation::test_energy(void)
{
    // Test void constructor
    test_try("Void constructor");
    try {
        GEnergy energy;
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test value constructors
    GEnergy erg(3.0, "erg");
    test_value(erg.erg(), 3.0, 1.0e-6, "GEnergy(erg) constructor");
    GEnergy keV(3.0, "keV");
    test_value(keV.keV(), 3.0, 1.0e-6, "GEnergy(keV) constructor");
    GEnergy MeV(3.0, "MeV");
    test_value(MeV.MeV(), 3.0, 1.0e-6, "GEnergy(MeV) constructor");
    GEnergy GeV(3.0, "GeV");
    test_value(GeV.GeV(), 3.0, 1.0e-6, "GEnergy(GeV) constructor");
    GEnergy TeV(3.0, "TeV");
    test_value(TeV.TeV(), 3.0, 1.0e-6, "GEnergy(TeV) constructor");
    GEnergy Angstrom(3.0, "Angstrom");
    test_value(Angstrom.Angstrom(), 3.0, 1.0e-6, "GEnergy(Angstrom) constructor");

    // Test energy value access methods
    test_value(MeV.erg(), 3.0/gammalib::erg2MeV, 1.0e-6, "erg() get method");
    test_value(MeV.keV(), 3.0e+3, 1.0e-6, "keV() get method");
    test_value(MeV.GeV(), 3.0e-3, 1.0e-6, "GeV() get method");
    test_value(MeV.TeV(), 3.0e-6, 1.0e-6, "TeV() get method");
    test_value(MeV.Angstrom(), gammalib::MeV2Angstrom/3.0, 1.0e-6, "Angstrom() get method");
    test_value(MeV("erg"), 3.0/gammalib::erg2MeV, 1.0e-6, "(erg) get operator");
    test_value(MeV("keV"), 3.0e+3, 1.0e-6, "(keV) get operator");
    test_value(MeV("MeV"), 3.0, 1.0e-6, "(MeV) get operator");
    test_value(MeV("GeV"), 3.0e-3, 1.0e-6, "(GeV) get operator");
    test_value(MeV("TeV"), 3.0e-6, 1.0e-6, "(TeV) get operator");
    test_value(MeV("Angstrom"), gammalib::MeV2Angstrom/3.0, 1.0e-6, "(Angstrom) get operator");
    test_value(MeV.log10erg(), std::log10(3.0/gammalib::erg2MeV), 1.0e-6, "log10erg() get method");
    test_value(MeV.log10keV(), std::log10(3.0e+3), 1.0e-6, "log10keV() get method");
    test_value(MeV.log10MeV(), std::log10(3.0), 1.0e-6, "log10MeV() get method");
    test_value(MeV.log10GeV(), std::log10(3.0e-3), 1.0e-6, "log10GeV() get method");
    test_value(MeV.log10TeV(), std::log10(3.0e-6), 1.0e-6, "log10TeV() get method");
    test_value(MeV.log10("erg"), std::log10(3.0/gammalib::erg2MeV), 1.0e-6, "log10(erg) get method");
    test_value(MeV.log10("keV"), std::log10(3.0e+3), 1.0e-6, "log10(keV) get method");
    test_value(MeV.log10("MeV"), std::log10(3.0), 1.0e-6, "log10(MeV) get method");
    test_value(MeV.log10("GeV"), std::log10(3.0e-3), 1.0e-6, "log10(GeV) get method");
    test_value(MeV.log10("TeV"), std::log10(3.0e-6), 1.0e-6, "log10(TeV) get method");

    // Test value set methods
    GEnergy energy;
    energy.erg(3.0);
    test_value(energy.erg(), 3.0, 1.0e-6, "erg() set method");
    energy.keV(3.0);
    test_value(energy.keV(), 3.0, 1.0e-6, "keV() set method");
    energy.MeV(3.0);
    test_value(energy.MeV(), 3.0, 1.0e-6, "MeV() set method");
    energy.GeV(3.0);
    test_value(energy.GeV(), 3.0, 1.0e-6, "GeV() set method");
    energy.TeV(3.0);
    test_value(energy.TeV(), 3.0, 1.0e-6, "TeV() set method");
    energy.Angstrom(3.0);
    test_value(energy.Angstrom(), 3.0, 1.0e-6, "Angstrom() set method");
    energy(3.0, "erg");
    test_value(energy.erg(), 3.0, 1.0e-6, "(erg) set operator");
    energy(3.0, "keV");
    test_value(energy.keV(), 3.0, 1.0e-6, "(keV) set method");
    energy(3.0, "MeV");
    test_value(energy.MeV(), 3.0, 1.0e-6, "(MeV) set method");
    energy(3.0, "GeV");
    test_value(energy.GeV(), 3.0, 1.0e-6, "(GeV) set method");
    energy(3.0, "TeV");
    test_value(energy.TeV(), 3.0, 1.0e-6, "(TeV) set method");
    energy(3.0, "Angstrom");
    test_value(energy.Angstrom(), 3.0, 1.0e-6, "(Angstrom) set method");
    energy.log10erg(std::log10(3.0));
    test_value(energy.erg(), 3.0, 1.0e-6, "log10erg() set method");
    energy.log10keV(std::log10(3.0));
    test_value(energy.keV(), 3.0, 1.0e-6, "log10keV() set method");
    energy.log10MeV(std::log10(3.0));
    test_value(energy.MeV(), 3.0, 1.0e-6, "log10MeV() set method");
    energy.log10GeV(std::log10(3.0));
    test_value(energy.GeV(), 3.0, 1.0e-6, "log10GeV() set method");
    energy.log10TeV(std::log10(3.0));
    test_value(energy.TeV(), 3.0, 1.0e-6, "log10TeV() set method");
    energy.log10(std::log10(3.0), "erg");
    test_value(energy.erg(), 3.0, 1.0e-6, "log10(erg) set method");
    energy.log10(std::log10(3.0), "keV");
    test_value(energy.keV(), 3.0, 1.0e-6, "log10(keV) set method");
    energy.log10(std::log10(3.0), "MeV");
    test_value(energy.MeV(), 3.0, 1.0e-6, "log10(MeV) set method");
    energy.log10(std::log10(3.0), "GeV");
    test_value(energy.GeV(), 3.0, 1.0e-6, "log10(GeV) set method");
    energy.log10(std::log10(3.0), "TeV");
    test_value(energy.TeV(), 3.0, 1.0e-6, "log10(TeV) set method");

    // Test operators
    double  va = 2.1;
    double  vb = 3.7;
    GEnergy a(va, "MeV");
    GEnergy b(vb, "MeV");
    GEnergy test = a + b;
    test_value(test.MeV(), va+vb, 1.0e-6, "operator+");
    test = a - b;
    test_value(test.MeV(), va-vb, 1.0e-6, "operator-");
    test = a * 2.9;
    test_value(test.MeV(), va*2.9, 1.0e-6, "operator*");
    test = 2.8 * a;
    test_value(test.MeV(), 2.8*va, 1.0e-6, "operator*");
    test = a / 3.7;
    test_value(test.MeV(), va/3.7, 1.0e-6, "operator*");
    test  = a;
    test += b;
    test_value(test.MeV(), va+vb, 1.0e-6, "operator+=");
    test  = a;
    test -= b;
    test_value(test.MeV(), va-vb, 1.0e-6, "operator-=");
    test_assert((a == a), "operator== (equal values)");
    test_assert(!(a == b), "operator== (non equal values)");
    test_assert(!(a != a), "operator!= (equal values)");
    test_assert((a != b), "operator!= (non equal values)");
    test_assert((a < b), "less than operator");
    test_assert(!(b < a), "less than operator");
    test_assert(!(a > b), "larger than operator");
    test_assert((b > a), "larger than operator");
    test_assert((a <= a), "less or equal operator");
    test_assert((a <= b), "less or equal operator");
    test_assert(!(b <= a), "less or equal operator");
    test_assert((a >= a), "larger or equal operator");
    test_assert(!(a >= b), "larger or equal operator");
    test_assert((b >= a), "larger or equal operator");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test GEnergies
 ***************************************************************************/
void TestGObservation::test_energies(void)
{
    // Test void constructor
    test_try("Void constructor");
    try {
        GEnergies energies;
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Manipulate GEnergies starting from an empty object
    GEnergies energies;
    test_value(energies.size(), 0, "GEnergies should have zero size.");
    test_assert(energies.is_empty(), "GEnergies should be empty.");

    // Add an energy
    energies.append(GEnergy());
    test_value(energies.size(), 1, "GEnergies should have 1 energy.");
    test_assert(!energies.is_empty(), "GEnergies should not be empty.");

    // Remove energy
    energies.remove(0);
    test_value(energies.size(), 0, "GEnergies should have zero size.");
    test_assert(energies.is_empty(), "GEnergies should be empty.");

    // Append two energies
    energies.append(GEnergy());
    energies.append(GEnergy());
    test_value(energies.size(), 2, "GEnergies should have 2 energies.");
    test_assert(!energies.is_empty(), "GEnergies should not be empty.");

    // Clear object
    energies.clear();
    test_value(energies.size(), 0, "GEnergies should have zero size.");
    test_assert(energies.is_empty(), "GEnergies should be empty.");

    // Insert two energies
    energies.insert(0, GEnergy());
    energies.insert(0, GEnergy());
    test_value(energies.size(), 2, "GEnergies should have 2 energies.");
    test_assert(!energies.is_empty(), "GEnergies should not be empty.");

    // Extend energies
    energies.extend(energies);
    test_value(energies.size(), 4, "GEnergies should have 4 energies.");
    test_assert(!energies.is_empty(), "GEnergies should not be empty.");

    // Create 4 energies
    energies.clear();
    for (int i = 0; i < 4; ++i) {
        energies.append(GEnergy(double(i), "MeV"));
    }
    for (int i = 0; i < 4; ++i) {
        test_value(energies[i].MeV(), double(i));
    }

    // Remove test file
    GFilename filename("test_energies.fits");
    filename.remove();

    // Save and reload energies
    test_try("Saving and loading");
    try {
        energies.save("test_energies.fits", true);
        energies.clear();
        energies.load("test_energies.fits");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }
    for (int i = 0; i < 4; ++i) {
        test_value(energies[i].MeV(), double(i));
    }

    // Test load constructor
    test_try("Load constructor");
    try {
        GEnergies energies2("test_energies.fits");
        for (int i = 0; i < 4; ++i) {
            test_value(energies[i].MeV(), double(i));
        }
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Save and reload energies in another extension
    test_try("Saving and loading in another extension");
    try {
        energies.save("test_energies.fits[NEW ENERGIES]", true);
        energies.clear();
        energies.load("test_energies.fits[NEW ENERGIES]");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }
    for (int i = 0; i < 4; ++i) {
        test_value(energies[i].MeV(), double(i));
    }

    // Check linear energies
    energies.clear();
    energies.set_lin(3, GEnergy(1.0, "MeV"), GEnergy(3.0, "MeV"));
    test_value(energies.size(), 3, "GEnergies should have 3 elements.");
    test_assert(!energies.is_empty(), "GEnergies should not be empty.");
    test_value(energies[0].MeV(), 1.0, 1.0e-10, "Energy 0 should be 1 MeV.");
    test_value(energies[1].MeV(), 2.0, 1.0e-10, "Energy 1 should be 2 MeV.");
    test_value(energies[2].MeV(), 3.0, 1.0e-10, "Energy 3 should be 3 MeV.");

    // Check logarithmic energies
    energies.clear();
    energies.set_log(3, GEnergy(1.0, "MeV"), GEnergy(100.0, "MeV"));
    test_value(energies.size(), 3, "GEnergies should have 3 elements.");
    test_assert(!energies.is_empty(), "GEnergies should not be empty.");
    test_value(energies[0].MeV(), 1.0, 1.0e-10, "Energy 0 should be 1 MeV.");
    test_value(energies[1].MeV(), 10.0, 1.0e-10, "Energy 1 should be 10 MeV.");
    test_value(energies[2].MeV(), 100.0, 1.0e-10, "Energy 2 should be 100 MeV.");

    // Check energy boundary set method
    GEbounds ebds1(2, GEnergy(1.0, "MeV"), GEnergy(3.0, "MeV"), false);
    energies.set(ebds1);
    test_value(energies.size(), 3, "GEbounds constructor (3 elements)");
    test_value(energies[0].MeV(), 1.0, 1.0e-10, "GEbounds constructor (3 elements)");
    test_value(energies[1].MeV(), 2.0, 1.0e-10, "GEbounds constructor (3 elements)");
    test_value(energies[2].MeV(), 3.0, 1.0e-10, "GEbounds constructor (3 elements)");

    // Check energy boundary set method
    GEbounds ebds2(1, GEnergy(1.0, "MeV"), GEnergy(1.0, "MeV"));
    energies.set(ebds2);
    test_value(energies.size(), 1, "GEbounds constructor (1 element)");
    test_value(energies[0].MeV(), 1.0, 1.0e-10, "GEbounds constructor (1 element)");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test GPhotons
 ***************************************************************************/
void TestGObservation::test_photons(void)
{
    // Test void constructor
    test_try("Void constructor");
    try {
        GPhotons photons;
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Manipulate GPhotons starting from an empty object
    GPhotons photons;
    test_value(photons.size(), 0, "GPhotons should have zero size.");
    test_assert(photons.is_empty(), "GPhotons should be empty.");

    // Add one photon
    photons.append(GPhoton());
    test_value(photons.size(), 1, "GPhotons should have 1 photon.");
    test_assert(!photons.is_empty(), "GPhotons should not be empty.");

    // Remove photon
    photons.remove(0);
    test_value(photons.size(), 0, "GPhotons should have zero size.");
    test_assert(photons.is_empty(), "GPhotons should be empty.");

    // Append two photons
    photons.append(GPhoton());
    photons.append(GPhoton());
    test_value(photons.size(), 2, "GPhotons should have 2 photons.");
    test_assert(!photons.is_empty(), "GPhotons should not be empty.");

    // Clear object
    photons.clear();
    test_value(photons.size(), 0, "GPhotons should have zero size.");
    test_assert(photons.is_empty(), "GPhotons should be empty.");

    // Insert two photons
    photons.insert(0, GPhoton());
    photons.insert(0, GPhoton());
    test_value(photons.size(), 2, "GPhotons should have 2 photons.");
    test_assert(!photons.is_empty(), "GPhotons should not be empty.");

    // Extend photons
    photons.extend(photons);
    test_value(photons.size(), 4, "GPhotons should have 4 photons.");
    test_assert(!photons.is_empty(), "GPhotons should not be empty.");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test GObservations::likelihood
 ***************************************************************************/
void TestGObservation::test_observations_optimizer(void)
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

    // Add observations to container
    for (int i = 0; i < 6; ++i) {

        // Set random Generator
        GRan ran(i);

        // Create an event list
        GEvents *events = model.generateList(rate, tmin, tmax, ran);

        // Create a observation
        GTestObservation observation;
        observation.id(gammalib::str(i));

        // Add events to the observation
        observation.events(*events);
        observation.ontime(tmax.secs()-tmin.secs());

        // Append observation to container
        obs.append(observation);

        // Free the event list
        delete events;
        
    } // endfor: added observations to container

    // Add the model to the observation
    obs.models(models);

    // Create a GLog for show the interations of optimizer.
    GLog log;

    // Create an optimizer.
    GOptimizerLM opt(&log);

    // Set number of stalls
    opt.max_stalls(50);

    // Optimize
    obs.optimize(opt);

    // Compute errors
    obs.errors(opt);

    // Get the result
    GModelPar result = (*(obs.models()[0]))[0];

    // Check if converged
    test_value(opt.status(), 0, "Check if converged",
                                "Optimizer did not converge");

    // Check if value is correct
    test_value(result.factor_value(), RATE, result.factor_error()*3.0);

    // Save covariance matrix as FITS file
    obs.function().save("test_observations_optimizer.fits");

    // Try loading covariance matrix as FITS file
    test_try("Load covariance matrix as FITS file");
    try {
        GFits file("test_observations_optimizer.fits");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Save covariance matrix as CSV file
    obs.function().save("test_observations_optimizer.csv");

    // Try loading covariance matrix as CSV file
    test_try("Load covariance matrix as CSV file");
    try {
        GCsv file("test_observations_optimizer.csv");
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

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

    // Append unbinned tests
    append(static_cast<pfunction>(&TestOpenMP::test_observations_optimizer_unbinned_1), "Test unbinned optimization (1 thread)");
    append(static_cast<pfunction>(&TestOpenMP::test_observations_optimizer_unbinned_10), "Test unbinned optimization (10 threads)");

    // Append binned tests
    append(static_cast<pfunction>(&TestOpenMP::test_observations_optimizer_binned_1), "Test binned optimization (1 thread)");
    append(static_cast<pfunction>(&TestOpenMP::test_observations_optimizer_binned_10), "Test binned optimisation (10 threads)");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone test suite
 *
 * @return Pointer to deep copy of test suite.
 ***************************************************************************/
TestOpenMP* TestOpenMP::clone(void) const
{
    // Clone test suite
    return new TestOpenMP(*this);
}


/***********************************************************************//**
 * @brief Test observations optimizer.
 *
 * @param[in] mode Testing mode.
 * 
 * This method supports two testing modes: 0 = unbinned and 1 = binned.
 ***************************************************************************/
void TestOpenMP::test_observations_optimizer(const int& mode)
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
        ob.events(*events);
        ob.ontime(tmax.secs()-tmin.secs());
        obs.append(ob);

        // Delete events pointer
        delete events;
        
    }

    // Add the model to the observation
    obs.models(models);

    // Create a GLog for show the interations of optimizer.
    GLog log;

    // Create an optimizer.
    GOptimizerLM opt(&log);

    // Set number of stalls
    opt.max_stalls(50);

    // Optimize
    obs.optimize(opt);

    // Compute errors
    obs.errors(opt);

    // Get the result
    GModelPar result = (*(obs.models()[0]))[0];

    // Check if converged
    test_assert(opt.status()==0, "Check if converged", 
                                 "Optimizer did not converge"); 

    // Check if value is correct
    test_value(result.factor_value(),RATE,result.factor_error()*3); 

    // Return
    return;
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
