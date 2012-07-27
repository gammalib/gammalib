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
#include <iostream>
#include <unistd.h>
#include "GCOMLib.hpp"
#include "GTools.hpp"

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
 * @brief Checks handling of IAQ response files
 *
 * This function checks the handling of IAQ response files. IAQ response
 * files are 2D images that show the instrument response as function of
 * geometrical (Phi_geo) and measured (Phi_bar) Compton scatter angle.
 ***************************************************************************/
void test_iaq_response(void)
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
 * @brief Checks handling of COMPTEL response files
 ***************************************************************************/
void test_response(void)
{
    // Dump header
    std::cout << "Test COMPTEL response: ";

    // Test IAQ
    test_iaq_response();

    // Dump final ok
    std::cout << " ok." << std::endl;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Checks handling of COMPTEL event bin
 ***************************************************************************/
void test_event_bin(void)
{
    // Dump header
    std::cout << "Test COMPTEL event bin: ";

    // Test COMPTEL event bin methods (one by one)
    try {
        // Event bin void constructor
        GCOMEventBin bin;
        std::cout << ".";

        // Copy constructor
        GCOMEventBin bin2(bin);
        std::cout << ".";

        // Assignment operator
        GCOMEventBin bin3 = bin;
        std::cout << ".";

        // clear method
        bin.clear();
        std::cout << ".";

        // clone method
        GCOMEventBin* bin4 = bin.clone();
        std::cout << ".";

        // size method
        if (bin.size() != 0) {
            std::cout << std::endl
                      << "TEST ERROR: Event bin size is expected to be zero"
                      << " but a size of "+str(bin.size())+" was found."
                      << std::endl;
            throw;
        };
        std::cout << ".";

        // dir method
        GCOMInstDir dir = bin.dir();
        std::cout << ".";

        // energy method
        GEnergy energy = bin.energy();
        std::cout << ".";

        // time method
        GTime time = bin.time();
        std::cout << ".";

        // counts method
        double counts = bin.counts();
        if (counts != 0) {
            std::cout << std::endl
                      << "TEST ERROR: Counts are expected to be zero"
                      << " but "+str(counts)+" counts were found."
                      << std::endl;
            throw;
        };
        std::cout << ".";

        // error method
        double error = bin.error();
        if (error > 1.0e-10) { // A small delta is added for the optimizer
            std::cout << std::endl
                      << "TEST ERROR: The counts error is expected to be zero"
                      << " but "+str(error)+" error counts were found."
                      << std::endl;
            throw;
        };
        std::cout << ".";

        // counts setting
        bin.counts(1.0);
        std::cout << ".";

        // print method
        std::string text = bin.print();
        std::cout << ".";

        // omega method
        double omega = bin.omega();
        if (omega != 0) {
            std::cout << std::endl
                      << "TEST ERROR: The solid angle is expected to be zero"
                      << " but an angle of "+str(omega)+" was found."
                      << std::endl;
            throw;
        };
        std::cout << ".";

        // ewidth method
        GEnergy ewidth = bin.ewidth();
        std::cout << ".";

        // ontime method
        double ontime = bin.ontime();
        if (ontime != 0) {
            std::cout << std::endl
                      << "TEST ERROR: The ontime is expected to be zero"
                      << " but an ontime of "+str(ontime)+" was found."
                      << std::endl;
            throw;
        };
        std::cout << ".";

    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to handle empty COMPTEL event bin."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }

    // Test COMPTEL event bin operations
    try {
        // Write and read back counts
        GCOMEventBin bin;
        bin.counts(3.3);
        std::cout << ".";
        double counts = bin.counts();
        if (counts != 3.3) {
            std::cout << std::endl
                      << "TEST ERROR: Counts are expected to be 3.3"
                      << " but "+str(counts)+" counts were found."
                      << std::endl;
            throw;
        };
        std::cout << ".";

    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to handle COMPTEL event bin."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }

    // Dump final ok
    std::cout << " ok." << std::endl;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Checks handling of COMPTEL event cube
 ***************************************************************************/
void test_event_cube(void)
{
    // Dump header
    std::cout << "Test COMPTEL event cube: ";

    // Test COMPTEL event cube methods (one by one)
    try {
        // Event cube void constructor
        GCOMEventCube cube;
        std::cout << ".";

        // Event cube load constructor
        GCOMEventCube cube2(com_dre);
        std::cout << ".";

        // Event cube copy constructor
        GCOMEventCube cube3(cube2);
        if (cube2.size() != cube3.size()) {
            std::cout << std::endl
                      << "TEST ERROR: Different cube size after using the"
                      << " copy constructor (before="+str(cube2.size())
                      << " after="+str(cube3.size())+")"
                      << std::endl;
            throw;
        };
        if (cube2.number() != cube3.number()) {
            std::cout << std::endl
                      << "TEST ERROR: Different number of events after using"
                      << " the copy constructor (before="+str(cube2.number())
                      << " after="+str(cube3.number())+")"
                      << std::endl;
            throw;
        };
        std::cout << ".";

        // Event cube assignment operator
        GCOMEventCube cube4 = cube2;
        if (cube2.size() != cube4.size()) {
            std::cout << std::endl
                      << "TEST ERROR: Different cube size after using the"
                      << " assignment operator (before="+str(cube2.size())
                      << " after="+str(cube4.size())+")"
                      << std::endl;
            throw;
        };
        if (cube2.number() != cube4.number()) {
            std::cout << std::endl
                      << "TEST ERROR: Different number of events after using"
                      << " the assignment operator (before="+str(cube2.number())
                      << " after="+str(cube4.number())+")"
                      << std::endl;
            throw;
        };
        std::cout << ".";
        
        // clear method
        cube4.clear();
        if (cube4.size() != 0) {
            std::cout << std::endl
                      << "TEST ERROR: Expected event cube with size 0 after"
                      << " clear but found size "+str(cube4.size())+"."
                      << std::endl;
            throw;
        };
        if (cube4.number() != 0) {
            std::cout << std::endl
                      << "TEST ERROR: Expected 0 events in cube after"
                      << " clear but found size "+str(cube4.number())+" events."
                      << std::endl;
            throw;
        };
        std::cout << ".";

        // clone method
        GCOMEventCube* cube5 = cube2.clone();
        if (cube2.size() != cube5->size()) {
            std::cout << std::endl
                      << "TEST ERROR: Different cube size after cloning"
                      << " (before="+str(cube2.size())
                      << " after="+str(cube5->size())+")"
                      << std::endl;
            throw;
        };
        if (cube2.number() != cube5->number()) {
            std::cout << std::endl
                      << "TEST ERROR: Different number of events after"
                      << " cloning (before="+str(cube2.number())
                      << " after="+str(cube5->number())+")"
                      << std::endl;
            throw;
        };
        std::cout << ".";

        // size method
        if (cube2.size() != 140600) {
            std::cout << std::endl
                      << "TEST ERROR: Expected cube dimension 140600, found "
                      << str(cube2.size())+"."
                      << std::endl;
            throw;
        };
        std::cout << ".";

        // dim method
        if (cube2.dim() != 3) {
            std::cout << std::endl
                      << "TEST ERROR: Expected cube dimension 3, found "
                      << str(cube2.dim())+"."
                      << std::endl;
            throw;
        };
        std::cout << ".";

        // naxis method
        if (cube2.naxis(0) != 76) {
            std::cout << std::endl
                      << "TEST ERROR: Expected Chi axis dimension 76, found "
                      << str(cube2.naxis(0))+"."
                      << std::endl;
            throw;
        };
        std::cout << ".";
        if (cube2.naxis(1) != 74) {
            std::cout << std::endl
                      << "TEST ERROR: Expected Psi axis dimension 74, found "
                      << str(cube2.naxis(1))+"."
                      << std::endl;
            throw;
        };
        std::cout << ".";
        if (cube2.naxis(2) != 25) {
            std::cout << std::endl
                      << "TEST ERROR: Expected Phi axis dimension 25, found "
                      << str(cube2.naxis(2))+"."
                      << std::endl;
            throw;
        };
        std::cout << ".";

        // nchi, npsi, nphi methods
        if (cube2.nchi() != 76) {
            std::cout << std::endl
                      << "TEST ERROR: Expected Chi axis dimension 76, found "
                      << str(cube2.nchi())+"."
                      << std::endl;
            throw;
        };
        std::cout << ".";
        if (cube2.npsi() != 74) {
            std::cout << std::endl
                      << "TEST ERROR: Expected Psi axis dimension 74, found "
                      << str(cube2.npsi())+"."
                      << std::endl;
            throw;
        };
        std::cout << ".";
        if (cube2.nphi() != 25) {
            std::cout << std::endl
                      << "TEST ERROR: Expected Phi axis dimension 25, found "
                      << str(cube2.nphi())+"."
                      << std::endl;
            throw;
        };
        std::cout << ".";

        // npix method
        int npix = cube2.nchi() * cube2.npsi();
        if (cube2.npix() != npix) {
            std::cout << std::endl
                      << "TEST ERROR: Expected "+str(npix)+" pixels in"
                      << " (Chi,Psi) plane, found "
                      << str(cube2.npix())+"."
                      << std::endl;
            throw;
        };

        // number method
        if (cube2.number() != 316141) {
            std::cout << std::endl
                      << "TEST ERROR: Expected 316141 events in cube, found "
                      << str(cube2.number())+"."
                      << std::endl;
            throw;
        };
        std::cout << ".";

        // event access
        double sum = 0.0;
        for (int i = 0; i < cube2.size(); ++i) {
            sum += cube2[i]->counts();
        }
        if (int(sum+0.5) != cube2.number()) {
            std::cout << std::endl
                      << "TEST ERROR: Expected "+str(cube2.number())
                      << " events in cube, found "
                      << str(int(sum+0.5))+" by summing over all elements."
                      << std::endl;
            throw;
        };
        std::cout << ".";

    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to handle COMPTEL event cube."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }

    // Dump final ok
    std::cout << " ok." << std::endl;

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
void test_binned_obs(void)
{
    // Write header
    std::cout << "Test binned COMPTEL observation handling: ";

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
    // Notify final test success
    std::cout << " ok." << std::endl;

    // Exit test
    return;
 
}


/***********************************************************************//**
 * @brief Main test function
 ***************************************************************************/
int main(void)
{
    // Dump header
    std::cout << std::endl;
    std::cout << "*********************************************" << std::endl;
    std::cout << "* COMPTEL instrument specific class testing *" << std::endl;
    std::cout << "*********************************************" << std::endl;

    // Set CALDB environment variable
    std::string caldb = "CALDB="+com_caldb;
    putenv((char*)caldb.c_str());

    // Execute tests not needing data
    test_event_bin();
    test_event_cube();
    test_response();
    test_binned_obs();

    // Return
    return 0;
}
