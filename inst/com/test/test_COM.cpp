/***************************************************************************
 *                      test_COM.cpp  -  test COM classes                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012 by Juergen Knoedlseder                              *
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
 * @file test_COM.cpp
 * @brief Testing of COM classes
 * @author <author>
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <stdlib.h>
#include <unistd.h>
#include "GTools.hpp"
#include "test_COM.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Globals ____________________________________________________________ */

/* __ Constants __________________________________________________________ */
const std::string datadir       = "../inst/com/test/data";
const std::string com_caldb     = "../inst/com/caldb";
const std::string com_iaq       = "u47569_iaq.fits";          // 1-3 MeV
const std::string com_dre       = datadir+"/m50439_dre.fits"; // 1-3 MeV
const std::string com_drb       = datadir+"/m34997_drg.fits";
const std::string com_drg       = datadir+"/m34997_drg.fits";
const std::string com_drx       = datadir+"/m32171_drx.fits";
const std::string com_obs       = datadir+"obs.xml";


/***********************************************************************//**
 * @brief Set COMPTEL response test methods
 ***************************************************************************/
void TestGCOMResponse::set(void)
{
    // Set test name
    name("GCOMResponse");

    // Append tests to test suite
    append(static_cast<pfunction>(&TestGCOMResponse::test_iaq_response), "Test IAQ response");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set COMPTEL observation test methods
 ***************************************************************************/
void TestGCOMObservation::set(void)
{
    // Set test name
    name("GCOMObservation");

    // Append tests to test suite
    append(static_cast<pfunction>(&TestGCOMObservation::test_binned_obs), "Test binned observation");
    append(static_cast<pfunction>(&TestGCOMObservation::test_event_bin), "Test event bin");
    append(static_cast<pfunction>(&TestGCOMObservation::test_event_cube), "Test event cube");

    // Return
    return;
}




/***********************************************************************//**
 * @brief Checks handling of IAQ response files
 *
 * This function checks the handling of IAQ response files. IAQ response
 * files are 2D images that show the instrument response as function of
 * geometrical (Phi_geo) and measured (Phi_bar) Compton scatter angle.
 ***************************************************************************/
void TestGCOMResponse::test_iaq_response(void)
{
/*
    // Test response loading
    try {
        // Construct observation from datasets
        GCOMResponse rsp(com_caldb, com_iaq);
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to construct IAQ response."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";
*/
    // Return
    return;
}


/***********************************************************************//**
 * @brief Checks handling of binned observations
 *
 * This function checks the handling of binned observations. Binned
 * observations are defined by
 * a DRE file containing the events,
 * a DRB file containing the background model,
 * a DRG file containing the geometrical probability of having an interaction
 *            in the second detector plane, and
 * a DRX file containing the sky exposure.
 *
 * Note that in first approximation we can use the DRG file as background
 * model.
 ***************************************************************************/
void TestGCOMObservation::test_binned_obs(void)
{
    // Declare observation and container
    GObservations obs;

/*
    // Dataset constructor
    try {
        // Construct observation from datasets
        GCOMObservation com(com_dre, com_drb, com_drg, com_drx);
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to construct observation from datasets."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // XML constructor
    try {
        // Construct observation from XML file
        GObservations obs(com_obs);
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to construct observation from XML file."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";
*/

    // Exit test
    return;
 
}


/***********************************************************************//**
 * @brief Checks handling of COMPTEL event bin
 ***************************************************************************/
void TestGCOMObservation::test_event_bin(void)
{
    // Test COMPTEL event bin methods (one by one)
    test_try("Test event bin methods");
    try {
        // Event bin void constructor
        GCOMEventBin bin;

        // Copy constructor
        GCOMEventBin bin2(bin);

        // Assignment operator
        GCOMEventBin bin3 = bin;

        // clear method
        bin.clear();

        // clone method
        GCOMEventBin* bin4 = bin.clone();

        // size method
        test_value(bin.size(), 0, "Test size() method.");

        // dir method
        GCOMInstDir dir = bin.dir();

        // energy method
        GEnergy energy = bin.energy();

        // time method
        GTime time = bin.time();

        // counts method
        test_value(bin.counts(), 0.0, 1.0e-10, "Test counts() method for zero counts.");

        // error method
        test_value(bin.error(), 0.0, 1.0e-10, "Test error() method for zero counts.");

        // counts setting
        bin.counts(1.0);
        test_value(bin.counts(), 1.0, 1.0e-10, "Test counts() method for one count.");
        test_value(bin.error(), 1.0, 1.0e-10, "Test error() method for one count.");

        // print method
        std::string text = bin.print();

        // omega method
        test_value(bin.omega(), 0.0, 1.0e-10, "Test omega() method.");

        // ewidth method
        GEnergy ewidth = bin.ewidth();

        // ontime method
        test_value(bin.ontime(), 0.0, 1.0e-10, "Test ontime() method.");

        // If we arrived here, signal success
        test_try_success();
    }
    catch (std::exception &e) {
        test_try_failure(e);
    }

    // Test COMPTEL event bin operations
    GCOMEventBin bin;
    bin.counts(3.3);
    test_value(bin.counts(), 3.3, 1.0e-10, "Test counts() method for 3.3 counts.");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Checks handling of COMPTEL event cube
 ***************************************************************************/
void TestGCOMObservation::test_event_cube(void)
{
    // Event cube void constructor
    GCOMEventCube cube;

    // Event cube load constructor
    GCOMEventCube cube2(com_dre);

    // Event cube copy constructor
    GCOMEventCube cube3(cube2);
    test_value(cube2.size(), cube3.size(), "Test copy constructor.",
               "Different cube size after using the"
               " copy constructor (before="+str(cube2.size())+
               " after="+str(cube3.size())+")");
    test_value(cube2.number(), cube3.number(), "Test copy constructor.",
               "Different number of events after using the"
               " copy constructor (before="+str(cube2.number())+
               " after="+str(cube3.number())+").");

    // Event cube assignment operator
    GCOMEventCube cube4 = cube2;
    test_value(cube2.size(), cube4.size(), "Test assignment operator.",
               "Different cube size after using the"
               " assignment operator (before="+str(cube2.size())+
               " after="+str(cube4.size())+")");
    test_value(cube2.number(), cube4.number(), "Test assignment operator.",
               "Different number of events after using the"
               " assignment operator (before="+str(cube2.number())+
               " after="+str(cube4.number())+").");
        
    // clear method
    cube4.clear();
    test_value(cube4.size(), 0, "Test clear method.",
               "Expected event cube with size 0 after"
               " clear but found size "+str(cube4.size())+".");
    test_value(cube4.number(), 0, "Test clear method.",
               "Expected 0 events in cube after"
               " clear but found "+str(cube4.size())+" events.");

    // clone method
    GCOMEventCube* cube5 = cube2.clone();
    test_value(cube2.size(), cube5->size(), "Test clone() method.",
               "Different cube size after cloning"
               " (before="+str(cube2.size())+
               " after="+str(cube5->size())+")");
    test_value(cube2.number(), cube5->number(), "Test clone() method.",
               "Different number of events after cloning"
               " (before="+str(cube2.number())+
               " after="+str(cube5->number())+").");

    // size method
    test_value(cube2.size(), 140600, "Test size() method.",
               "Expected cube dimension 140600, found "+
               str(cube2.size())+".");

    // dim method
    test_value(cube2.dim(), 3, "Test dim() method.",
               "Expected 3 cube dimensions, found "+
               str(cube2.dim())+".");

    // naxis method
    test_value(cube2.naxis(0), 76, "Test naxis(0) method.",
               "Expected Chi axis dimension 76, found "+
               str(cube2.naxis(0))+".");
    test_value(cube2.naxis(1), 74, "Test naxis(1) method.",
               "Expected Chi axis dimension 74, found "+
               str(cube2.naxis(1))+".");
    test_value(cube2.naxis(2), 25, "Test naxis(2) method.",
               "Expected Chi axis dimension 25, found "+
               str(cube2.naxis(2))+".");

    // nchi, npsi, nphi methods
    test_value(cube2.nchi(), 76, "Test nchi() method.",
               "Expected Chi axis dimension 76, found "+
               str(cube2.nchi())+".");
    test_value(cube2.npsi(), 74, "Test npsi() method.",
               "Expected Chi axis dimension 74, found "+
               str(cube2.npsi())+".");
    test_value(cube2.nphi(), 25, "Test nphi() method.",
               "Expected Chi axis dimension 25, found "+
               str(cube2.nphi())+".");

    // npix method
    int npix = cube2.nchi() * cube2.npsi();
    test_value(cube2.npix(), npix, "Test npix() method.",
               "Expected "+str(npix)+" pixels in (Chi,Psi) plane, found "+
               str(cube2.npix())+".");

    // number method
    test_value(cube2.number(), 316141, "Test number() method.",
               "Expected 316141 events in cube, found "+
               str(cube2.number())+".");

    // event access
    double sum = 0.0;
    for (int i = 0; i < cube2.size(); ++i) {
        sum += cube2[i]->counts();
    }
    test_value(cube2.number(), int(sum+0.5), "Test event access.",
               "Expected "+str(cube2.number())+
               " events in cube, found "+
               str(int(sum+0.5))+" by summing over all elements.");

    // Return
    return;
}


/***************************************************************************
 * @brief Main entry point for test executable
 ***************************************************************************/
int main(void)
{
    // Allocate test suit container
    GTestSuites testsuites("COMPTEL instrument specific class testing");

    // Set CALDB environment variable
    std::string caldb = "CALDB="+com_caldb;
    putenv((char*)caldb.c_str());

    // Initially assume that we pass all tests
    bool success = true;

    // Create test suites and append them to the container
    TestGCOMResponse    rsp;
    TestGCOMObservation obs;
    testsuites.append(rsp);
    testsuites.append(obs);

    // Run the testsuites
    success = testsuites.run();

    // Save test report
    testsuites.save("reports/GCOM.xml");

    // Return success status
    return (success ? 0 : 1);
}
