/***************************************************************************
 *              test_GObservation.cpp - Test observation module            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2022 by Jean-Baptiste Cayrou                        *
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
//#include <cstdlib>                 // std::getenv
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

/* __ Namespaces _________________________________________________________ */

/* __ Globals ____________________________________________________________ */

/* __ Constants __________________________________________________________ */
const std::string datadir        = gammalib::getenv("TEST_DATA");
const std::string ephem_integral = datadir + "/ephem_crab_integral.fits";
const std::string ephem_tempo2   = datadir + "/ephem_crab_tempo2.par";
const std::string ephem_fermi    = datadir + "/ephem_db_fermi.fits";
const std::string ephem_psrtime  = datadir + "/ephem_xte.psrtime";

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
    append(static_cast<pfunction>(&TestGObservation::test_pulsar),
           "Test GPulsar class");
    append(static_cast<pfunction>(&TestGObservation::test_pulsar_ephemeris),
           "Test GPulsarEphemeris class");
    append(static_cast<pfunction>(&TestGObservation::test_ephemerides),
           "Test GEphemerides class");
    append(static_cast<pfunction>(&TestGObservation::test_phases),
           "Test GPhases class");
    append(static_cast<pfunction>(&TestGObservation::test_photons),
           "Test GPhotons class");
    append(static_cast<pfunction>(&TestGObservation::test_response_cache),
           "Test GResponseCache class");
    append(static_cast<pfunction>(&TestGObservation::test_response_vector_cache),
           "Test GResponseVectorCache class");
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

    // Remove fully enclosed energy interval
    ebds.clear();
    ebds.append(GEnergy(1.0, "MeV"), GEnergy(10.0, "MeV"));
    ebds.append(GEnergy(10.0, "MeV"), GEnergy(100.0, "MeV"));
    ebds.remove(GEnergy(2.0, "MeV"), GEnergy(3.0, "MeV"));
    test_value(ebds.size(), 3, "GEbounds should have 3 element.");
    test_value(ebds.emin(0).MeV(), 1.0, 1.0e-10, "First energy minimum should be 1.");
    test_value(ebds.emax(0).MeV(), 2.0, 1.0e-10, "First energy maximum should be 2.");
    test_value(ebds.emin(1).MeV(), 3.0, 1.0e-10, "Second energy minimum should be 3.");
    test_value(ebds.emax(1).MeV(), 10.0, 1.0e-10, "Second energy minimum should be 10.");

    // Remove matching energy interval
    ebds.clear();
    ebds.append(GEnergy(1.0, "MeV"), GEnergy(10.0, "MeV"));
    ebds.append(GEnergy(10.0, "MeV"), GEnergy(100.0, "MeV"));
    ebds.remove(GEnergy(1.0, "MeV"), GEnergy(10.0, "MeV"));
    test_value(ebds.size(), 1, "GEbounds should have 1 element.");
    test_value(ebds.emin(0).MeV(), 10.0, 1.0e-10, "First energy minimum should be 10.");
    test_value(ebds.emax(0).MeV(), 100.0, 1.0e-10, "First energy maximum should be 100.");

    // Remove overlapping energy interval
    ebds.clear();
    ebds.append(GEnergy(1.0, "MeV"), GEnergy(10.0, "MeV"));
    ebds.append(GEnergy(10.0, "MeV"), GEnergy(100.0, "MeV"));
    ebds.remove(GEnergy(8.0, "MeV"), GEnergy(12.0, "MeV"));
    test_value(ebds.size(), 2, "GEbounds should have 2 element.");
    test_value(ebds.emin(0).MeV(), 1.0, 1.0e-10, "First energy minimum should be 1.");
    test_value(ebds.emax(0).MeV(), 8.0, 1.0e-10, "First energy maximum should be 8.");
    test_value(ebds.emin(1).MeV(), 12.0, 1.0e-10, "Second energy minimum should be 12.");
    test_value(ebds.emax(1).MeV(), 100.0, 1.0e-10, "Second energy minimum should be 100.");

    // Remove partially overlapping energy interval
    ebds.clear();
    ebds.append(GEnergy(1.0, "MeV"), GEnergy(10.0, "MeV"));
    ebds.append(GEnergy(10.0, "MeV"), GEnergy(100.0, "MeV"));
    ebds.remove(GEnergy(80.0, "MeV"), GEnergy(200.0, "MeV"));
    test_value(ebds.size(), 2, "GEbounds should have 2 element.");
    test_value(ebds.emin(0).MeV(), 1.0, 1.0e-10, "First energy minimum should be 1.");
    test_value(ebds.emax(0).MeV(), 10.0, 1.0e-10, "First energy maximum should be 10.");
    test_value(ebds.emin(1).MeV(), 10.0, 1.0e-10, "Second energy minimum should be 10.");
    test_value(ebds.emax(1).MeV(), 80.0, 1.0e-10, "Second energy minimum should be 80.");

    // Remove partially overlapping energy interval
    ebds.clear();
    ebds.append(GEnergy(1.0, "MeV"), GEnergy(10.0, "MeV"));
    ebds.append(GEnergy(10.0, "MeV"), GEnergy(100.0, "MeV"));
    ebds.remove(GEnergy(0.1, "MeV"), GEnergy(2.0, "MeV"));
    test_value(ebds.size(), 2, "GEbounds should have 2 element.");
    test_value(ebds.emin(0).MeV(), 2.0, 1.0e-10, "First energy minimum should be 2.");
    test_value(ebds.emax(0).MeV(), 10.0, 1.0e-10, "First energy maximum should be 10.");
    test_value(ebds.emin(1).MeV(), 10.0, 1.0e-10, "Second energy minimum should be 10.");
    test_value(ebds.emax(1).MeV(), 100.0, 1.0e-10, "Second energy minimum should be 100.");

    // Remove non overlapping energy interval
    ebds.clear();
    ebds.append(GEnergy(1.0, "MeV"), GEnergy(10.0, "MeV"));
    ebds.append(GEnergy(10.0, "MeV"), GEnergy(100.0, "MeV"));
    ebds.remove(GEnergy(100.0, "MeV"), GEnergy(1000.0, "MeV"));
    test_value(ebds.size(), 2, "GEbounds should have 2 element.");
    test_value(ebds.emin(0).MeV(), 1.0, 1.0e-10, "First energy minimum should be 1.");
    test_value(ebds.emax(0).MeV(), 10.0, 1.0e-10, "First energy maximum should be 10.");
    test_value(ebds.emin(1).MeV(), 10.0, 1.0e-10, "Second energy minimum should be 10.");
    test_value(ebds.emax(1).MeV(), 100.0, 1.0e-10, "Second energy minimum should be 100.");

    // Check linear boundaries
    ebds.clear();
    ebds.set(2, GEnergy(1.0, "MeV"), GEnergy(3.0, "MeV"), "LIN");
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
    ebds.set(2, GEnergy(1.0, "MeV"), GEnergy(100.0, "MeV"), "LOG");
    test_value(ebds.size(), 2, "GEbounds should have 2 elements.");
    test_assert(!ebds.is_empty(), "GEbounds should not be empty.");
    test_value(ebds.emin(0).MeV(), 1.0, 1.0e-10, "Bin 0 minimum energy should be 1.");
    test_value(ebds.emin(1).MeV(), 10.0, 1.0e-10, "Bin 1 minimum energy should be 10.");
    test_value(ebds.emax(0).MeV(), 10.0, 1.0e-10, "Bin 0 maximum energy should be 10.");
    test_value(ebds.emax(1).MeV(), 100.0, 1.0e-10, "Bin 1 maximum energy should be 100.");
    test_value(ebds.emin().MeV(), 1.0, 1.0e-10, "Minimum energy should be 1.");
    test_value(ebds.emax().MeV(), 100.0, 1.0e-10, "Maximum energy should be 100.");

    // Check power-law boundaries
    ebds.clear();
    ebds.set(2, GEnergy(1.0, "MeV"), GEnergy(100.0, "MeV"), "POW", 1.0);
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
    ebds.set(1, GEnergy(1.0, "MeV"), GEnergy(10.0, "MeV"), "LOG");
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
        GEbounds bad(10, GEnergy(100.0, "MeV"), GEnergy(10.0, "MeV"), "LIN");
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
        GEbounds bad(10, GEnergy(100.0, "MeV"), GEnergy(10.0, "MeV"), "LOG");
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
        GEbounds bad(10, GEnergy(0.0, "MeV"), GEnergy(10.0, "MeV"), "LOG");
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
        GEbounds bad(10, GEnergy(10.0, "MeV"), GEnergy(0.0, "MeV"), "LOG");
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
    GEnergies energies1(3, GEnergy(1.0, "MeV"), GEnergy(10.0, "MeV"), "LIN");
    GEbounds ebds1(energies1);
    test_value(ebds1.size(), 2, "GEnergies constructor (3 elements)");
    test_value(ebds1.emin().MeV(),  1.0, 1.0e-10, "GEnergies constructor (3 elements)");
    test_value(ebds1.emax(0).MeV(), 5.5, 1.0e-10, "GEnergies constructor (3 elements)");
    test_value(ebds1.emax().MeV(), 10.0, 1.0e-10, "GEnergies constructor (3 elements)");

    // Check energies constructor method
    GEnergies energies2(1, GEnergy(1.0, "MeV"), GEnergy(1.0, "MeV"), "LIN");
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

    // Check equality operator
    test_assert(ebds1 == ebds1, "Equality operator");

    // Check inequality operator
    test_assert(ebds1 != ebds2, "Inequality operator");

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
    test_value(gti.ontime(), 0.0, 1.0e-10, "Check ontime of empty object");
    test_value(gti.telapse(), 0.0, 1.0e-10, "Check elapsed time of empty object");

    // Add empty interval
    gti.append(GTime(1.0), GTime(1.0));
    test_value(gti.size(), 1, "GGti should have 1 interval.");
    test_assert(!gti.is_empty(), "GGti should not be empty.");
    test_value(gti.tstart().secs(), 1.0, 1.0e-7, "Start time should be 1.");
    test_value(gti.tstop().secs(), 1.0, 1.0e-7, "Stop time should be 1.");
    test_value(gti.ontime(), 0.0, 1.0e-10, "Check ontime of empty interval");
    test_value(gti.telapse(), 0.0, 1.0e-10, "Check elapsed time of empty interval");

    // Add one interval
    gti.append(GTime(1.0), GTime(10.0));
    test_value(gti.size(), 2, "GGti should have 2 intervals.");
    test_assert(!gti.is_empty(), "GGti should not be empty.");
    test_value(gti.tstart().secs(), 1.0, 1.0e-7, "Start time should be 1.");
    test_value(gti.tstop().secs(), 10.0, 1.0e-7, "Stop time should be 10.");
    test_value(gti.ontime(), 9.0, 1.0e-10, "Check ontime of one interval");
    test_value(gti.telapse(), 9.0, 1.0e-10, "Check elapsed time of one interval");

    // Remove interval
    gti.remove(0);
    gti.remove(0);
    test_value(gti.size(), 0, "GGti should have zero size.");
    test_assert(gti.is_empty(), "GGti should be empty.");
    test_value(gti.tstart().secs(), 0.0, 1.0e-7, "Start time should be 0.");
    test_value(gti.tstop().secs(), 0.0, 1.0e-7, "Stop time should be 0.");
    test_value(gti.ontime(), 0.0, 1.0e-10, "Check ontime after removal");
    test_value(gti.telapse(), 0.0, 1.0e-10, "Check elapsed time after removal");

    // Append two overlapping intervals
    gti.append(GTime(1.0), GTime(100.0));
    gti.append(GTime(10.0), GTime(1000.0));
    test_value(gti.size(), 2, "GGti should have 2 intervals.");
    test_assert(!gti.is_empty(), "GGti should not be empty.");
    test_value(gti.tstart().secs(), 1.0, 1.0e-7, "Start time should be 1.");
    test_value(gti.tstop().secs(), 1000.0, 1.0e-7, "Stop time should be 1000.");
    test_value(gti.ontime(), 1089.0, 1.0e-10, "Check ontime of two overlapping intervals");
    test_value(gti.telapse(), 999.0, 1.0e-10, "Check elapsed time of two overlapping intervals");

    // Clear object
    gti.clear();
    test_value(gti.size(), 0, "GGti should have zero size.");
    test_assert(gti.is_empty(), "GGti should be empty.");
    test_value(gti.tstart().secs(), 0.0, 1.0e-7, "Start time should be 0.");
    test_value(gti.tstop().secs(), 0.0, 1.0e-7, "Stop time should be 0.");
    test_value(gti.ontime(), 0.0, 1.0e-10, "Check ontime after clear");
    test_value(gti.telapse(), 0.0, 1.0e-10, "Check elapsed time after clear");

    // Append two overlapping intervals in inverse order
    gti.clear();
    gti.append(GTime(10.0), GTime(1000.0));
    gti.append(GTime(1.0), GTime(100.0));
    test_value(gti.size(), 2, "GGti should have 2 intervals.");
    test_assert(!gti.is_empty(), "GGti should not be empty.");
    test_value(gti.tstart().secs(), 1.0, 1.0e-7, "Start time should be 1.");
    test_value(gti.tstop().secs(), 1000.0, 1.0e-7, "Stop time should be 1000.");
    test_value(gti.ontime(), 1089.0, 1.0e-10, "Check ontime of two overlapping intervals");
    test_value(gti.telapse(), 999.0, 1.0e-10, "Check elapsed time of two overlapping intervals");

    // Insert three non-overlapping intervals
    gti.clear();
    gti.insert(GTime(1.0), GTime(2.0));
    gti.insert(GTime(10.0), GTime(20.0));
    gti.insert(GTime(100.0), GTime(200.0));
    test_value(gti.size(), 3, "GGti should have 3 intervals.");
    test_assert(!gti.is_empty(), "GGti should not be empty.");
    test_value(gti.tstart().secs(), 1.0, 1.0e-10, "Start time should be 1.");
    test_value(gti.tstop().secs(), 200.0, 1.0e-10, "Stop time should be 200.");
    test_value(gti.ontime(), 111.0, 1.0e-10, "Check ontime of 3 non-overlapping intervals");
    test_value(gti.telapse(), 199.0, 1.0e-10, "Check elapsed time of 3 non-overlapping intervals");

    // Check merging
    gti.merge();
    test_value(gti.size(), 3, "Check size after merging");
    test_assert(!gti.is_empty(), "Check that GTI is not empty after merging");
    test_value(gti.tstart().secs(), 1.0, 1.0e-10, "Check start time after merging");
    test_value(gti.tstop().secs(), 200.0, 1.0e-10, "Check stop time after merging");
    test_value(gti.ontime(), 111.0, 1.0e-10, "Check ontime after merging");
    test_value(gti.telapse(), 199.0, 1.0e-10, "Check elapsed time after merging");

    // Insert three non-overlapping intervals in reverse order
    gti.clear();
    gti.insert(GTime(100.0), GTime(200.0));
    gti.insert(GTime(10.0), GTime(20.0));
    gti.insert(GTime(1.0), GTime(2.0));
    test_value(gti.size(), 3, "GGti should have 3 intervals.");
    test_assert(!gti.is_empty(), "GGti should not be empty.");
    test_value(gti.tstart().secs(), 1.0, 1.0e-10, "Start time should be 1.");
    test_value(gti.tstop().secs(), 200.0, 1.0e-10, "Stop time should be 200.");
    test_value(gti.ontime(), 111.0, 1.0e-10, "Check ontime of 3 non-overlapping intervals");
    test_value(gti.telapse(), 199.0, 1.0e-10, "Check elapsed time of 3 non-overlapping intervals");

    // Check merging
    gti.merge();
    test_value(gti.size(), 3, "Check size after merging");
    test_assert(!gti.is_empty(), "Check that GTI is not empty after merging");
    test_value(gti.tstart().secs(), 1.0, 1.0e-10, "Check start time after merging");
    test_value(gti.tstop().secs(), 200.0, 1.0e-10, "Check stop time after merging");
    test_value(gti.ontime(), 111.0, 1.0e-10, "Check ontime after merging");
    test_value(gti.telapse(), 199.0, 1.0e-10, "Check elapsed time after merging");

    // Insert three non-overlapping intervals in mixed order
    gti.clear();
    gti.insert(GTime(100.0), GTime(200.0));
    gti.insert(GTime(1.0), GTime(2.0));
    gti.insert(GTime(10.0), GTime(20.0));
    test_value(gti.size(), 3, "GGti should have 3 intervals.");
    test_assert(!gti.is_empty(), "GGti should not be empty.");
    test_value(gti.tstart().secs(), 1.0, 1.0e-10, "Start time should be 1.");
    test_value(gti.tstop().secs(), 200.0, 1.0e-10, "Stop time should be 200.");
    test_value(gti.ontime(), 111.0, 1.0e-10, "Check ontime of 3 non-overlapping intervals");
    test_value(gti.telapse(), 199.0, 1.0e-10, "Check elapsed time of 3 non-overlapping intervals");

    // Check merging
    gti.merge();
    test_value(gti.size(), 3, "Check size after merging");
    test_assert(!gti.is_empty(), "Check that GTI is not empty after merging");
    test_value(gti.tstart().secs(), 1.0, 1.0e-10, "Check start time after merging");
    test_value(gti.tstop().secs(), 200.0, 1.0e-10, "Check stop time after merging");
    test_value(gti.ontime(), 111.0, 1.0e-10, "Check ontime after merging");
    test_value(gti.telapse(), 199.0, 1.0e-10, "Check elapsed time after merging");

    // Insert two overlapping intervals
    gti.clear();
    gti.insert(GTime(1.0), GTime(100.0));
    gti.insert(GTime(10.0), GTime(1000.0));
    test_value(gti.size(), 2, "GGti should have 2 intervals.");
    test_assert(!gti.is_empty(), "GGti should not be empty.");
    test_value(gti.tstart().secs(), 1.0, 1.0e-10, "Start time should be 1.");
    test_value(gti.tstop().secs(), 1000.0, 1.0e-10, "Stop time should be 1000.");
    test_value(gti.ontime(), 1089.0, 1.0e-10, "Check ontime of 2 overlapping intervals");
    test_value(gti.telapse(), 999.0, 1.0e-10, "Check elapsed time of 2 overlapping intervals");

    // Check merging
    gti.merge();
    test_value(gti.size(), 1, "Check size after merging");
    test_assert(!gti.is_empty(), "Check that GTI is not empty after merging");
    test_value(gti.tstart().secs(), 1.0, 1.0e-10, "Check start time after merging");
    test_value(gti.tstop().secs(), 1000.0, 1.0e-10, "Check stop time after merging");
    test_value(gti.ontime(), 999.0, 1.0e-10, "Check ontime after merging");
    test_value(gti.telapse(), 999.0, 1.0e-10, "Check elapsed time after merging");

    // Insert two overlapping intervals in reverse order
    gti.clear();
    gti.insert(GTime(10.0), GTime(1000.0));
    gti.insert(GTime(1.0), GTime(100.0));
    test_value(gti.size(), 2, "GGti should have 2 intervals.");
    test_assert(!gti.is_empty(), "GGti should not be empty.");
    test_value(gti.tstart().secs(), 1.0, 1.0e-7, "Start time should be 1.");
    test_value(gti.tstop().secs(), 1000.0, 1.0e-7, "Stop time should be 1000.");
    test_value(gti.ontime(), 1089.0, 1.0e-10, "Check ontime of 2 overlapping intervals");
    test_value(gti.telapse(), 999.0, 1.0e-10, "Check elapsed time of 2 overlapping intervals");

    // Merge two overlapping intervals
    gti.clear();
    gti.merge(GTime(1.0), GTime(100.0));
    gti.merge(GTime(10.0), GTime(1000.0));
    test_value(gti.size(), 1, "GGti should have 1 interval.");
    test_assert(!gti.is_empty(), "GGti should not be empty.");
    test_value(gti.tstart().secs(), 1.0, 1.0e-7, "Start time should be 1.");
    test_value(gti.tstop().secs(), 1000.0, 1.0e-7, "Stop time should be 1000.");
    test_value(gti.ontime(), 999.0, 1.0e-10, "Check ontime of merged intervals");
    test_value(gti.telapse(), 999.0, 1.0e-10, "Check elapsed time of merged intervals");

    // Merge two overlapping intervals in inverse order
    gti.clear();
    gti.merge(GTime(10.0), GTime(1000.0));
    gti.merge(GTime(1.0), GTime(100.0));
    test_value(gti.size(), 1, "GGti should have 1 interval.");
    test_assert(!gti.is_empty(), "GGti should not be empty.");
    test_value(gti.tstart().secs(), 1.0, 1.0e-7, "Start time should be 1.");
    test_value(gti.tstop().secs(), 1000.0, 1.0e-7, "Stop time should be 1000.");
    test_value(gti.ontime(), 999.0, 1.0e-10, "Check ontime of merged intervals");
    test_value(gti.telapse(), 999.0, 1.0e-10, "Check elapsed time of merged intervals");

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

    // Check saving in a different extension
    test3.save("test_gti.fits[GOOD TIME INTERVALS]", true);
    GGti load2("test_gti.fits[GOOD TIME INTERVALS]");
    test_value(load2.size(), 1, "GGti should have 1 interval.");
    test_value(load2.tstart().convert(GTimeReference(51544.5, "s")), 0.0, 1.0e-7,
               "Start time should be 0.");
    test_value(load2.tstop().convert(GTimeReference(51544.5, "s")), 100.0, 1.0e-7,
               "Stop time should be 100.");

    // Check contains() method
    GGti test6;
    test6.append(GTime(2.0), GTime(3.0));
    test6.append(GTime(5.0), GTime(7.0));
    test6.append(GTime(10.0), GTime(20.0));
    test_value(test6.size(), 3, "Check GGti::contains(): size=3");
    test_value(test6.tstart().secs(), 2.0, 1.0e-7, "Check GGti::contains(): start=2");
    test_value(test6.tstop().secs(), 20.0, 1.0e-7, "Check GGti::contains(): stop=20");
    test_assert(!test6.contains(GTime(1.0)), "Check GGti::contains(): 1 is not contained");
    test_assert(test6.contains(GTime(15.0)), "Check GGti::contains(): 15 is contained");
    test_assert(test6.contains(GTime(20.0)), "Check GGti::contains(): 20 is contained");
    test_assert(test6.contains(GTime(2.5)), "Check GGti::contains(): 2.5 is contained");
    test_assert(test6.contains(GTime(20.0)), "Check GGti::contains(): 20 is still contained");
    test_assert(test6.contains(GTime(6.0)), "Check GGti::contains(): 6 is contained");
    test_assert(!test6.contains(GTime(8.0)), "Check GGti::contains(): 8 is not contained");
    test_assert(!test6.contains(GTime(4.0)), "Check GGti::contains(): 4 is not contained");

    // Check overlap() method
    GGti test7;
    test7.append(GTime(5.0), GTime(6.0));
    test7.append(GTime(9.0), GTime(10.0));
    test_value(test7.overlap(GTime(0.0), GTime(5.0)), 0.0, 1.0e-7,
               "Check GGti::overlap(): interval before GTIs");
    test_value(test7.overlap(GTime(10.0), GTime(15.0)), 0.0, 1.0e-7,
               "Check GGti::overlap(): interval after GTIs");
    test_value(test7.overlap(GTime(4.0), GTime(5.5)), 0.5, 1.0e-7,
               "Check GGti::overlap(): interval ends in GTI");
    test_value(test7.overlap(GTime(9.5), GTime(11.0)), 0.5, 1.0e-7,
               "Check GGti::overlap(): interval starts in GTI");
    test_value(test7.overlap(GTime(4.0), GTime(7.0)), 1.0, 1.0e-7,
               "Check GGti::overlap(): interval starts before and ends after GTI");
    test_value(test7.overlap(GTime(5.1), GTime(5.9)), 0.8, 1.0e-7,
               "Check GGti::overlap(): interval within GTI");
    test_value(test7.overlap(GTime(5.0), GTime(10.0)), 2.0, 1.0e-7,
               "Check GGti::overlap(): interval comprises exactly GTIs");
    test_value(test7.overlap(GTime(4.0), GTime(11.0)), 2.0, 1.0e-7,
               "Check GGti::overlap(): interval comprises GTIs");
    test_value(test7.overlap(GTime(5.5), GTime(9.5)), 1.0, 1.0e-7,
               "Check GGti::overlap(): interval partly overlaps GTIs");

    // Check reduce(GGti) method
    GGti test8;
    test8.append(GTime(1.0), GTime(3.0));
    test8.append(GTime(5.0), GTime(6.0));
    GGti test9 = test8;
    GGti test10;
    test10.append(GTime(1.0), GTime(6.0));
    test9.reduce(test10);
    test_value(test9.size(), 2, "Check GGti::reduce(): size=2");
    test_value(test9.ontime(), 3.0, "Check GGti::reduce(): ontime=3");
    test_value(test9.tstart(0).secs(), 1.0, 1.0e-7, "Check GGti::reduce(): tstart[0]=1");
    test_value(test9.tstop(0).secs(),  3.0, 1.0e-7, "Check GGti::reduce(): tstop[0]=3");
    test_value(test9.tstart(1).secs(), 5.0, 1.0e-7, "Check GGti::reduce(): tstart[1]=5");
    test_value(test9.tstop(1).secs(),  6.0, 1.0e-7, "Check GGti::reduce(): tstop[1]=6");
    test9 = test8;
    test10.clear();
    test10.append(GTime(0.0), GTime(10.0));
    test9.reduce(test10);
    test_value(test9.size(), 2, "Check GGti::reduce(): size=2");
    test_value(test9.ontime(), 3.0, "Check GGti::reduce(): ontime=3");
    test_value(test9.tstart(0).secs(), 1.0, 1.0e-7, "Check GGti::reduce(): tstart[0]=1");
    test_value(test9.tstop(0).secs(),  3.0, 1.0e-7, "Check GGti::reduce(): tstop[0]=3");
    test_value(test9.tstart(1).secs(), 5.0, 1.0e-7, "Check GGti::reduce(): tstart[1]=5");
    test_value(test9.tstop(1).secs(),  6.0, 1.0e-7, "Check GGti::reduce(): tstop[1]=6");
    test9 = test8;
    test10.clear();
    test10.append(GTime(0.0), GTime(4.0));
    test9.reduce(test10);
    test_value(test9.size(), 1, "Check GGti::reduce(): size=1");
    test_value(test9.ontime(), 2.0, "Check GGti::reduce(): ontime=2");
    test_value(test9.tstart(0).secs(), 1.0, 1.0e-7, "Check GGti::reduce(): tstart[0]=1");
    test_value(test9.tstop(0).secs(),  3.0, 1.0e-7, "Check GGti::reduce(): tstop[0]=3");
    test9 = test8;
    test10.clear();
    test10.append(GTime(1.5), GTime(2.0));
    test10.append(GTime(2.2), GTime(2.5));
    test10.append(GTime(5.1), GTime(5.2));
    test9.reduce(test10);
    test_value(test9.size(), 3, "Check GGti::reduce(): size=3");
    test_value(test9.ontime(), 0.9, "Check GGti::reduce(): ontime=0.9");
    test_value(test9.tstart(0).secs(), 1.5, 1.0e-7, "Check GGti::reduce(): tstart[0]=1.5");
    test_value(test9.tstop(0).secs(),  2.0, 1.0e-7, "Check GGti::reduce(): tstop[0]=2");
    test_value(test9.tstart(1).secs(), 2.2, 1.0e-7, "Check GGti::reduce(): tstart[1]=2.2");
    test_value(test9.tstop(1).secs(),  2.5, 1.0e-7, "Check GGti::reduce(): tstop[1]=2.5");
    test_value(test9.tstart(2).secs(), 5.1, 1.0e-7, "Check GGti::reduce(): tstart[1]=5.1");
    test_value(test9.tstop(2).secs(),  5.2, 1.0e-7, "Check GGti::reduce(): tstop[1]=5.2");

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

    // Test time constructor (string)
    test_try("Time constructor (string)");
    try {
        GTime time("1800.01");
        test_try_success();
        test_value(time.secs(), 1800.01);
        test_value(time.days(), 1800.01/86400.0);
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
    test_value(time.utc(), "2010-01-02T10:17:37");
    test_value(time.utc(2), "2010-01-02T10:17:36.79");

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

    // Test time string set method
    time.set("2005-10-08T14:30:25");
    test_value(time.utc(), "2005-10-08T14:30:25");
    time.set("123456.789");
    test_value(time.secs(), 123456.789);
    time.set("123456.789 (TT)");
    test_value(time.secs(), 123456.789);
    time.set("123424.605(TAI)");
    test_value(time.secs(), 123456.789);
    time.set("123390.60487(UTC)");
    test_value(time.secs(), 123456.789);
    time.set("MJD55198.42966648");
    test_value(time.mjd(), 55198.42966648);
    time.set("MJD 55198.42966648(TT)");
    test_value(time.mjd(), 55198.42966648);
    time.set("MJD 55198.42929398 (TAI)");
    test_value(time.mjd(), 55198.42966648);
    time.set("MJD 55198.42890046 (UTC)");
    test_value(time.mjd(), 55198.42966648);
    time.set("JD 2455198.92966648");
    test_value(time.jd(), 2455198.92966648);
    time.set("JD2455198.92966648(TT)");
    test_value(time.jd(), 2455198.92966648);
    time.set("JD 2455198.92929398 (TAI)");
    test_value(time.jd(), 2455198.92966648);
    time.set("JD  2455198.92890046 (UTC)");
    test_value(time.jd(), 2455198.92966648);

    // Test time string set method using reference system
    time.set("123456.789", GTimeReference(mjd_ref, "secs", "TT", "LOCAL"));
    test_value(time.days(), 1.4288979);
    time.set("1.4288979", GTimeReference(mjd_ref, "days", "TT", "LOCAL"));
    test_value(time.days(), 1.4288979);
    time.set("1.4288979 (TT)", GTimeReference(mjd_ref, "days", "UTC", "LOCAL"));
    test_value(time.days(), 1.4288979);
    time.set("1.4285254 (TAI)", GTimeReference(mjd_ref, "days", "TT", "LOCAL"));
    test_value(time.days(), 1.4288979);
    time.set("1.42813188 (UTC)", GTimeReference(mjd_ref, "days", "TT", "LOCAL"));
    test_value(time.days(), 1.4288979);

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

    // Test leap_seconds method
    test_value(time.leap_seconds(), 36.0, "Test leap_seconds() method");

    // Test utc2tt method
    test_value(time.utc2tt(), 68.184, "Test utc2tt() method");

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
    energies.set(3, GEnergy(1.0, "MeV"), GEnergy(3.0, "MeV"), "LIN");
    test_value(energies.size(), 3, "GEnergies should have 3 elements.");
    test_assert(!energies.is_empty(), "GEnergies should not be empty.");
    test_value(energies[0].MeV(), 1.0, 1.0e-10, "Energy 0 should be 1 MeV.");
    test_value(energies[1].MeV(), 2.0, 1.0e-10, "Energy 1 should be 2 MeV.");
    test_value(energies[2].MeV(), 3.0, 1.0e-10, "Energy 3 should be 3 MeV.");

    // Check logarithmic energies
    energies.clear();
    energies.set(3, GEnergy(1.0, "MeV"), GEnergy(100.0, "MeV"), "LOG");
    test_value(energies.size(), 3, "GEnergies should have 3 elements.");
    test_assert(!energies.is_empty(), "GEnergies should not be empty.");
    test_value(energies[0].MeV(), 1.0, 1.0e-10, "Energy 0 should be 1 MeV.");
    test_value(energies[1].MeV(), 10.0, 1.0e-10, "Energy 1 should be 10 MeV.");
    test_value(energies[2].MeV(), 100.0, 1.0e-10, "Energy 2 should be 100 MeV.");

    // Check power-law energies
    energies.clear();
    energies.set(3, GEnergy(1.0, "MeV"), GEnergy(100.0, "MeV"), "POW", 1.0);
    test_value(energies.size(), 3, "GEnergies should have 3 elements.");
    test_assert(!energies.is_empty(), "GEnergies should not be empty.");
    test_value(energies[0].MeV(), 1.0, 1.0e-10, "Energy 0 should be 1 MeV.");
    test_value(energies[1].MeV(), 10.0, 1.0e-10, "Energy 1 should be 10 MeV.");
    test_value(energies[2].MeV(), 100.0, 1.0e-10, "Energy 2 should be 100 MeV.");

    // Check energy boundary set method
    GEbounds ebds1(2, GEnergy(1.0, "MeV"), GEnergy(3.0, "MeV"), "LIN");
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
 * @brief Test GPulsar class
 ***************************************************************************/
void TestGObservation::test_pulsar(void)
{
    // Test void constructor
    GPulsar pulsar1;
    test_assert(pulsar1.is_empty(), "Test for empty void GPulsar object");
    test_value(pulsar1.size(), 0, "Test for no ephemerides in void GPulsar object");

    // Test filename constructor using INTEGRAL file
    GPulsar pulsar2(ephem_integral);
    test_assert(!pulsar2.is_empty(),
         "Test for non-empty GPulsar object build from INTEGRAL file");
    test_value(pulsar2.size(), 5,
         "Test for 5 ephemerides in GPulsar object build from INTEGRAL file");
    test_value(pulsar2.validity().tstart().mjd(), 52670.0,
         "Test validity start for ephemerides in GPulsar object build from INTEGRAL file");
    test_value(pulsar2.validity().tstop().mjd(), 53290.0,
         "Test validity stop for ephemerides in GPulsar object build from INTEGRAL file");
    const GPulsarEphemeris& ephemeris2 = pulsar2[0];
    test_value(ephemeris2.name(), "Crab",
         "Test name of pulsar in INTEGRAL ephemeris file");
    test_value(ephemeris2.tstart().mjd(), 52670.0,
         "Test MJD start of ephemeris in INTEGRAL file");
    test_value(ephemeris2.tstop().mjd(), 52699.0,
         "Test MJD stop of ephemeris in INTEGRAL file");
    test_value(ephemeris2.timesys(), "TT",
         "Test timesystem of ephemeris in INTEGRAL file");
    test_value(ephemeris2.t0().mjd(), 52686.0,
         "Test epoch of ephemeris in INTEGRAL file");
    test_value(ephemeris2.phase(), 0.29227057582842,
         "Test phase of ephemeris in INTEGRAL file");
    test_value(ephemeris2.f0(), 29.8092382284677,
         "Test F0 of ephemeris in INTEGRAL file");
    test_value(ephemeris2.f1(), -3.736594e-10,
         "Test F1 of ephemeris in INTEGRAL file");
    test_value(ephemeris2.f2(), 2.9131e-20,
         "Test F2 of ephemeris in INTEGRAL file");

    // Test filename constructor using tempo2 file
    GPulsar pulsar3(ephem_tempo2);
    test_assert(!pulsar3.is_empty(),
         "Test for non-empty GPulsar object build from tempo2 file");
    test_value(pulsar3.size(), 1,
         "Test for 1 ephemeris in GPulsar object build from tempo2 file");
    test_value(pulsar3.validity().tstart().mjd(), 54588.676894515279855,
         "Test validity start for ephemerides in GPulsar object build from tempo2 file");
    test_value(pulsar3.validity().tstop().mjd(), 55257.704889561078744,
         "Test validity stop for ephemerides in GPulsar object build from tempo2 file");
    const GPulsarEphemeris& ephemeris3 = pulsar3[0];
    test_value(ephemeris3.name(), "PSR J0534+2200",
         "Test name of pulsar in tempo2 ephemeris file");
    test_value(ephemeris3.tstart().mjd(), 54588.676894515279855,
         "Test MJD start of ephemeris in tempo2 file");
    test_value(ephemeris3.tstop().mjd(), 55257.704889561078744,
         "Test MJD stop of ephemeris in tempo2 file");
    test_value(ephemeris3.timesys(), "TT",
         "Test timesystem of ephemeris in tempo2 file");
    test_value(ephemeris3.t0().mjd(), 54673.454536766465928,
         "Test epoch of ephemeris in tempo2 file");
    test_value(ephemeris3.phase(), 0.0,
         "Test phase of ephemeris in tempo2 file");
    test_value(ephemeris3.f0(), 29.745209371530433057,
         "Test F0 of ephemeris in tempo2 file");
    test_value(ephemeris3.f1(), -3.7194575283370609575e-10,
         "Test F1 of ephemeris in tempo2 file");
    test_value(ephemeris3.f2(), 1.1593944148649215896e-20,
         "Test F2 of ephemeris in tempo2 file");

    // Test filename constructor using Fermi file
    GPulsar pulsar4(ephem_fermi, "PSR J1028-5820");
    test_assert(!pulsar4.is_empty(),
         "Test for non-empty GPulsar object build from Fermi file");
    test_value(pulsar4.size(), 1,
         "Test for 1 ephemeris in GPulsar object build from Fermi file");
    test_value(pulsar4.validity().tstart().mjd(), 54563.0,
         "Test validity start for ephemerides in GPulsar object build from Fermi file");
    test_value(pulsar4.validity().tstop().mjd(), 54786.0,
         "Test validity stop for ephemerides in GPulsar object build from Fermi file");
    const GPulsarEphemeris& ephemeris4 = pulsar4[0];
    test_value(ephemeris4.name(), "PSR J1028-5820",
         "Test name of pulsar in Fermi ephemeris file");
    test_value(ephemeris4.tstart().mjd(), 54563.0,
         "Test MJD start of ephemeris in Fermi file");
    test_value(ephemeris4.tstop().mjd(), 54786.0,
         "Test MJD stop of ephemeris in Fermi file");
    test_value(ephemeris4.timesys(), "TT",
         "Test timesystem of ephemeris in Fermi file");
    test_value(ephemeris4.t0().mjd(), 54564.0,
         "Test epoch of ephemeris in Fermi file");
    test_value(ephemeris4.phase(), 0.955045061768033,
         "Test phase of ephemeris in Fermi file");
    test_value(ephemeris4.f0(), 10.9405324785477,
         "Test F0 of ephemeris in Fermi file");
    test_value(ephemeris4.f1(), -1.92828681288825e-12,
         "Test F1 of ephemeris in Fermi file");
    test_value(ephemeris4.f2(), 0.0,
         "Test F2 of ephemeris in Fermi file");

    // Test filename constructor using psrtime file
    GPulsar pulsar5(ephem_psrtime, "0531+21");
    test_assert(!pulsar5.is_empty(),
         "Test for non-empty GPulsar object build from psrtime file");
    test_value(pulsar5.size(), 56,
         "Test for 56 ephemerides in GPulsar object build from psrtime file");
    test_value(pulsar5.validity().tstart().mjd(), 48282.0,
         "Test validity start for ephemerides in GPulsar object build from psrtime file");
    test_value(pulsar5.validity().tstop().mjd(), 51122.0,
         "Test validity stop for ephemerides in GPulsar object build from psrtime file");
    const GPulsarEphemeris& ephemeris5 = pulsar5[0];
    test_value(ephemeris5.name(), "PSR B0531+21",
         "Test name of pulsar in psrtime ephemeris file");
    test_value(ephemeris5.tstart().mjd(), 48282.0,
         "Test MJD start of ephemeris in psrtime file");
    test_value(ephemeris5.tstop().mjd(), 48316.0,
         "Test MJD stop of ephemeris in psrtime file");
    test_value(ephemeris5.timesys(), "UTC",
         "Test timesystem of ephemeris in psrtime file");
    test_value(ephemeris5.t0().mjd(), 48299.0,
         "Test epoch of ephemeris in psrtime file");
    test_value(ephemeris5.phase(), -0.845097897454252,
         "Test phase of ephemeris in psrtime file");
    test_value(ephemeris5.f0(), 29.9516010684378,
         "Test F0 of ephemeris in psrtime file");
    test_value(ephemeris5.f1(), -3.77726e-10,
         "Test F1 of ephemeris in psrtime file");
    test_value(ephemeris5.f2(), 0.0,
         "Test F2 of ephemeris in psrtime file");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test GPulsarEphemeris class
 ***************************************************************************/
void TestGObservation::test_pulsar_ephemeris(void)
{
    // Test void constructor
    GPulsarEphemeris ephemeris1;
    test_value(ephemeris1.name(), "", "Test name() of void GPulsarEphemeris object");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test GEphemerides class
 ***************************************************************************/
void TestGObservation::test_ephemerides(void)
{
    // Test void constructor
    GEphemerides ephemerides1;
    test_assert(ephemerides1.is_empty(), "Test for empty void GEphemerides object");
    test_value(ephemerides1.size(), 0, "Test for no ephemerides in void GEphemerides object");

    // Test automatic ephemerides fetching procedure
    GEphemerides ephemerides2;
    test_assert(ephemerides2.is_empty(),
         "Test for empty GEphemerides object before call to geo2sbb");
    test_value(ephemerides2.size(), 0,
         "Test for no ephemerides in GEphemerides object before call to geo2sbb");
    test_value(ephemerides2.name(), "",
         "Test name for no ephemerides in GEphemerides object before call to geo2sbb");
    GSkyDir dir;
    GTime   time;
    GVector vector(3);
    double geo2ssb = ephemerides2.geo2ssb(dir, time);
    test_assert(!ephemerides2.is_empty(),
         "Test for non-empty GEphemerides object after call to geo2sbb");
    test_value(ephemerides2.size(), 32894,
         "Test for 32894 ephemerides in GEphemerides object after call to geo2sbb");
    test_value(ephemerides2.name(), "DE200",
         "Test name for ephemerides in GEphemerides object after call to geo2sbb");
    test_value(geo2ssb, -89.7124092753901, "Test for geo2sbb value");
    geo2ssb = ephemerides2.geo2ssb(dir, time, vector);
    test_assert(!ephemerides2.is_empty(),
         "Test for non-empty GEphemerides object after call to geo2sbb");
    test_value(ephemerides2.size(), 32894,
         "Test for 32894 ephemerides in GEphemerides object after call to geo2sbb");
    test_value(ephemerides2.name(), "DE200",
         "Test name for ephemerides in GEphemerides object after call to geo2sbb");
    test_value(geo2ssb, -89.7124092753901, "Test for geo2sbb value");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test GPhases class
 ***************************************************************************/
void TestGObservation::test_phases(void)
{
    // Test void constructor
    test_try("Void constructor");
    try {
        GPhases phases;
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Manipulate GPhases starting from an empty object
    GPhases phases;
    test_value(phases.size(), 0);
    test_assert(phases.is_empty(), "Test for empty GPhases object");

    // Add one phase interval
    phases.append(0.3, 0.5);
    test_value(phases.size(), 1);
    test_value(phases.pmin(0), 0.3);
    test_value(phases.pmax(0), 0.5);
    test_assert(!phases.is_empty(), "Test for non-empty GPhases object");

    // Remove phase interval
    phases.remove(0);
    test_value(phases.size(), 0);
    test_assert(phases.is_empty(), "Test for empty GPhases object");

    // Append two phase intervals
    phases.append(0.3, 0.5);
    phases.append(0.4, 0.7);
    test_value(phases.pmin(0), 0.3);
    test_value(phases.pmax(0), 0.5);
    test_value(phases.pmin(1), 0.4);
    test_value(phases.pmax(1), 0.7);
    test_value(phases.size(), 2);
    test_assert(!phases.is_empty(), "Test for non-empty GPhases object");

    // Copy phase intervals
    GPhases phases2(phases);
    test_value(phases2.pmin(0), 0.3);
    test_value(phases2.pmax(0), 0.5);
    test_value(phases2.pmin(1), 0.4);
    test_value(phases2.pmax(1), 0.7);
    test_value(phases2.size(), 2);
    test_assert(!phases2.is_empty(), "Test for non-empty GPhases object");

    // Clear object
    phases.clear();
    test_value(phases.size(), 0);
    test_assert(phases.is_empty(), "Test for empty GPhases object");

    // Extend photons
    phases.append(0.3, 0.5);
    phases.append(0.4, 0.7);
    phases.extend(phases);
    test_value(phases.size(), 4);
    test_assert(!phases.is_empty(), "Test for non-empty GPhases object");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test GResponseCache
 ***************************************************************************/
void TestGObservation::test_response_cache(void)
{
    // Allocate cache
    GResponseCache cache;

    // Set a few energies
    GEnergy      eng10(1.0, "MeV");
    GEnergy      eng12(1.2, "MeV");
    GEnergy      eng21(2.1, "MeV");
    GEnergy      eng37(3.7, "MeV");
    GTestInstDir dir1;
    GTestInstDir dir2;
    GTestInstDir dir3;
    dir1.hash(1);
    dir2.hash(2);
    dir3.hash(3);

    // Test empty cache
    test_assert(cache.is_empty(), "Test that void cache is empty");
    test_value(cache.size(), 0, "Test void cache size");

    // Add one element and check
    cache.set("Crab", eng10, eng12, 1010.0);
    test_assert(!cache.is_empty(), "Test that filled cache is not empty");
    test_value(cache.size(), 1, "Test filled cache size");
    double value = 0.0;
    bool flag = cache.contains("Crab", eng10, eng12, &value);
    test_assert(flag, "Test that cache contains Crab/eng10/eng12 element.");
    test_value(value, 1010.0, 1.0e-6, "Test cache value");

    // Add another element and check
    cache.set("Crab", dir1, eng10, eng21, 1021.0);
    test_value(cache.size(), 2, "Test filled cache size");
    flag = cache.contains("Crab", dir1, eng10, eng21, &value);
    test_assert(flag, "Test that cache contains Crab/eng10/eng21 element.");
    test_value(value, 1021.0, 1.0e-6, "Test cache value");

    // Add a few more elements
    cache.set("Crab", dir1, eng10, eng21, 1021.0); // Was already in, check update
    cache.set("Crab", dir3, eng10, eng37, 1037.0);
    cache.set("Crab", dir1, eng12, eng10, 1210.0);
    cache.set("Crab", dir2, eng12, eng12, 1212.0);
    cache.set("Crab", dir2, eng12, eng21, 1221.0);
    cache.set("Crab", dir2, eng12, eng37, 1237.0);
    cache.set("Vela", dir2, eng10, eng37, 1037.0);

    // Test cache now
    test_value(cache.size(), 8, "Test filled cache size");
    flag = cache.contains("Crab", dir1, eng10, eng21, &value);
    test_assert(flag, "Test that cache contains Crab/eng10/eng21 element.");
    test_value(value, 1021.0, 1.0e-6, "Test cache value");
    flag = cache.contains("Crab", dir2, eng12, eng12, &value);
    test_assert(flag, "Test that cache contains Crab/eng12/eng12 element.");
    test_value(value, 1212.0, 1.0e-6, "Test cache value");
    flag = cache.contains("Vela", dir2, eng10, eng37, &value);
    test_assert(flag, "Test that cache contains Crab/eng10/eng37 element.");
    test_value(value, 1037.0, 1.0e-6, "Test cache value");

    // Now copy the cache and test the copy
    GResponseCache copy_cache = cache;
    test_value(copy_cache.size(), 8, "Test filled copied cache size");
    flag = copy_cache.contains("Crab", dir2, eng12, eng12, &value);
    test_assert(flag, "Test that copied cache contains Crab/eng12/eng12 element.");
    test_value(value, 1212.0, 1.0e-6, "Test copied cache value");
    flag = copy_cache.contains("Vela", dir2, eng10, eng37, &value);
    test_assert(flag, "Test that copied cache contains Crab/eng10/eng37 element.");
    test_value(value, 1037.0, 1.0e-6, "Test copied cache value");

    // Now remove the Crab values and test the remaining cache
    cache.remove("Crab");
    test_value(cache.size(), 1, "Test cache size after Crab removal");
    flag = cache.contains("Vela", dir2, eng10, eng37, &value);
    test_assert(flag, "Test that cache contains Vela/eng10/eng37 element.");
    test_value(value, 1037.0, 1.0e-6, "Test cache value");

    // Now clear the cache and test that it is empty
    cache.clear();
    test_assert(cache.is_empty(), "Test that emptied cache is empty");
    test_value(cache.size(), 0, "Test emptied cache size");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test GResponseCache
 ***************************************************************************/
void TestGObservation::test_response_vector_cache(void)
{
    // Allocate cache
    GResponseVectorCache cache;

    // Set a few energies
    GVector vector1(10);
    vector1[1] = 1.0;
    vector1[2] = 2.0;
    GVector vector2(10);
    vector2[5] = 5.0;
    vector2[7] = 7.0;
    GVector vector3(12);
    vector3[1]  = 1.0;
    vector3[11] = 11.0;

    // Test empty cache
    test_assert(cache.is_empty(), "Test that void cache is empty");
    test_value(cache.size(), 0, "Test void cache size");

    // Add one vector and check
    cache.set("Crab", vector1);
    test_assert(!cache.is_empty(), "Test that filled cache is not empty");
    test_value(cache.size(), 1, "Test filled cache size");
    GVector vector4(10);
    bool flag = cache.contains("Crab", &vector4);
    test_assert(flag, "Test that cache contains Crab vector.");
    test_assert(vector4 == vector1, "Test cache values for Crab.");

    // Add another vector and check
    cache.set("Vela", vector2);
    test_assert(!cache.is_empty(), "Test that filled cache is not empty");
    test_value(cache.size(), 2, "Test filled cache size");
    GVector vector5(10);
    flag = cache.contains("Vela", &vector5);
    test_assert(flag, "Test that cache contains Vela vector.");
    test_assert(vector5 == vector2, "Test cache values for Vela.");

    // Replace Crab values (same dimension)
    cache.set("Crab", vector2);
    test_assert(!cache.is_empty(), "Test that filled cache is not empty");
    test_value(cache.size(), 2, "Test filled cache size");
    GVector vector6(10);
    flag = cache.contains("Crab", &vector6);
    test_assert(flag, "Test that cache contains Crab vector.");
    test_assert(vector6 == vector2, "Test cache values for replaced Crab.");

    // Replace Crab values (different dimension)
    cache.set("Crab", vector3);
    test_assert(!cache.is_empty(), "Test that filled cache is not empty");
    test_value(cache.size(), 2, "Test filled cache size");
    GVector vector7(12);
    flag = cache.contains("Crab", &vector7);
    test_assert(flag, "Test that cache contains Crab vector.");
    test_assert(vector7 == vector3, "Test cache values for extended Crab.");

    // Remove Crab cache
    cache.remove("Crab");
    test_value(cache.size(), 1, "Test filled cache size");
    flag = cache.contains("Crab", &vector7);
    test_assert(!flag, "Test that cache does not contain Crab vector.");
    flag = cache.contains("Vela", &vector5);
    test_assert(flag, "Test that cache contains Vela vector.");
    test_assert(vector5 == vector2, "Test cache values for Vela.");

    // Now clear the cache and test that it is empty
    cache.clear();
    test_assert(cache.is_empty(), "Test that emptied cache is empty");
    test_value(cache.size(), 0, "Test emptied cache size");

    // Add Crab and Vela to response vector cache and test saving and loading
    cache.set("Crab", vector1);
    cache.set("Vela", vector2);
    cache.save("test_response_vector_cache.fits", true);
    cache.clear();
    cache.load("test_response_vector_cache.fits");
    test_assert(!cache.is_empty(), "Test that loaded cache is not empty");
    test_value(cache.size(), 2, "Test loaded cache size");
    flag = cache.contains("Crab", &vector4);
    test_assert(flag, "Test that loaded cache contains Crab vector.");
    test_assert(vector4 == vector1, "Test loaded cache values for Crab.");
    flag = cache.contains("Vela", &vector5);
    test_assert(flag, "Test that loaded cache contains Vela vector.");
    test_assert(vector5 == vector2, "Test loaded cache values for Vela.");

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

    // Check results
    test_value(obs.nobserved(),   140415,   "Check number of observed events");
    test_value(obs.npred(),       140415.0, "Check number of predicted events");
    test_value(obs.npred("Test"), 140415.0, "Check number of predicted events for model");

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

    // Create and append test suites
    TestGObservation obs;
    testsuites.append(obs);

    // Create and append OpenMP test suite
    #ifdef _OPENMP
    TestOpenMP openmp;
    testsuites.append(openmp);
    #endif

    // Run the testsuites
    bool success = testsuites.run();

    // Save test report
    testsuites.save("reports/GObservation.xml");

    // Return success status
    return (success ? 0 : 1);
}
