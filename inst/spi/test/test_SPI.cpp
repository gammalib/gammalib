/***************************************************************************
 *                 test_SPI.cpp - Test INTEGRAL/SPI classes                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2020 by Juergen Knoedlseder                              *
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
 * @file test_SPI.cpp
 * @brief INTEGRAL/SPI test class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdlib>     // getenv
#include "test_SPI.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Globals ____________________________________________________________ */

/* __ Constants __________________________________________________________ */
const std::string datadir   = std::getenv("TEST_SPI_DATA");
const std::string spi_caldb = datadir + "/../../caldb";
const std::string og_dol    = datadir+"/obs/og_spi.fits";


/***********************************************************************//**
 * @brief Set INTEGRAL/SPI test methods
 ***************************************************************************/
void TestGSPI::set(void)
{
    // Set test name
    name("GSPI");

    // Append tests to test suite
    append(static_cast<pfunction>(&TestGSPI::test_instdir),
           "Test GSPIInstDir");
    append(static_cast<pfunction>(&TestGSPI::test_eventbin),
           "Test GSPIEventBin");
    append(static_cast<pfunction>(&TestGSPI::test_eventcube),
           "Test GSPIEventCube");
    append(static_cast<pfunction>(&TestGSPI::test_modeldataspace),
           "Test GSPIModelDataSpace");
    append(static_cast<pfunction>(&TestGSPI::test_obs),
           "Test GSPIObservation");
    append(static_cast<pfunction>(&TestGSPI::test_response),
           "Test GSPIResponse");
    append(static_cast<pfunction>(&TestGSPI::test_tools),
           "Test GSPITools");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone test suite
 *
 * @return Pointer to deep copy of test suite.
 ***************************************************************************/
TestGSPI* TestGSPI::clone(void) const
{
    // Clone test suite
    return new TestGSPI(*this);
}


/***********************************************************************//**
 * @brief Test GSPIInstDir class
 ***************************************************************************/
void TestGSPI::test_instdir(void)
{
    // Test content of void class
    GSPIInstDir dir1;
    test_value(dir1.classname(), "GSPIInstDir", "Void instance classname");
    test_value(dir1.dir().ra_deg(), 0.0, "Void instance Right Ascension");
    test_value(dir1.dir().dec_deg(), 0.0, "Void instance Declination");
    test_value(dir1.detid(), 0, "Void instance detector identifier");

    // Test content of filled class
    GSkyDir dir;
    dir.radec_deg(83.6331, 22.0145);
    GSPIInstDir dir2(dir, 13);
    test_value(dir2.dir().ra_deg(), 83.6331, "Filled instance Right Ascension");
    test_value(dir2.dir().dec_deg(), 22.0145, "Filled instance Declination");
    test_value(dir2.detid(), 13, "Filled instance detector identifier");

    // Test content of copied class
    GSPIInstDir dir3 = dir2;
    test_value(dir3.dir().ra_deg(), 83.6331, "Copied instance Right Ascension");
    test_value(dir3.dir().dec_deg(), 22.0145, "Copied instance Declination");
    test_value(dir3.detid(), 13, "Copied instance detector identifier");

    // Test content of copy-constructed class
    GSPIInstDir dir4(dir3);
    test_value(dir4.dir().ra_deg(), 83.6331, "Copied instance Right Ascension");
    test_value(dir4.dir().dec_deg(), 22.0145, "Copied instance Declination");
    test_value(dir4.detid(), 13, "Copied instance detector identifier");

    // Test setters and getters
    dir.radec_deg(299.5903, 35.2016);
    dir4.dir(dir);
    dir4.detid(1);
    test_value(dir4.dir().ra_deg(), 299.5903, "Get Right Ascension");
    test_value(dir4.dir().dec_deg(), 35.2016, "Get Declination");
    test_value(dir4.detid(), 1, "Get detector identifier");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test GSPIEventBin class
 ***************************************************************************/
void TestGSPI::test_eventbin(void)
{
    // Test content of void bin
    GSPIEventBin bin1;
    test_value(bin1.classname(), "GSPIEventBin", "Void instance classname");
    test_value(bin1.dir().dir().ra_deg(), 0.0, "Void instance Right Ascension");
    test_value(bin1.dir().dir().dec_deg(), 0.0, "Void instance Declination");
    test_value(bin1.dir().detid(), 0, "Void instance detector identifier");
    test_value(bin1.energy().keV(), 0.0, "Void instance energy");
    test_value(bin1.time().secs(), 0.0, "Void instance time");
    test_value(bin1.counts(), 0.0, "Void instance counts");
    test_value(bin1.error(), 0.0, "Void instance error");
    test_value(bin1.ontime(), 0.0, "Void instance ontime");
    test_value(bin1.livetime(), 0.0, "Void instance livetime");
    test_value(bin1.index(), -1, "Void instance index");
    test_value(bin1.ipt(), -1, "Void instance pointing index");
    test_value(bin1.idir(), -1, "Void instance direction index");
    test_value(bin1.iebin(), -1, "Void instance energy bin index");

    // Test counts setter and getter
    bin1.counts(41.0);
    test_value(bin1.counts(), 41.0, "Counts getter");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test GSPIEventCube class
 ***************************************************************************/
void TestGSPI::test_eventcube(void)
{
    // Test content of void cube
    GSPIEventCube cube1;
    test_value(cube1.classname(), "GSPIEventCube", "Void instance classname");
    test_value(cube1.size(), 0, "Void instance size");
    test_value(cube1.dim(), 3, "Void instance dimensions");
    test_value(cube1.naxis(0), 0, "Void instance pointings");
    test_value(cube1.naxis(1), 0, "Void instance detectors");
    test_value(cube1.naxis(2), 0, "Void instance energy bins");
    test_value(cube1.number(), 0, "Void instance number of events");
    test_value(cube1.ontime(), 0.0, "Void instance ontime");
    test_value(cube1.livetime(), 0.0, "Void instance livetime");
    test_value(cube1.models(), 0, "Void instance models");

    // Test loading of Observation Group
    cube1.load(og_dol);
    test_value(cube1.size(), 16720, "Loaded instance size");
    test_value(cube1.dim(), 3, "Loaded instance dimensions");
    test_value(cube1.naxis(0), 88, "Loaded instance pointings");
    test_value(cube1.naxis(1), 19, "Loaded instance detectors");
    test_value(cube1.naxis(2), 10, "Loaded instance energy bins");
    test_value(cube1.number(), 101269457, "Loaded instance number of events");
    test_value(cube1.ontime(), 193966.8178673, "Loaded instance ontime");
    test_value(cube1.livetime(), 170657.5371606, "Loaded instance livetime");
    test_value(cube1.models(), 1, "Loaded instance models");
    test_value(cube1.model_counts(0), 101269456.964783, "Loaded instance model counts");
    test_value(cube1.ptid(7), "00440010.000000", "Loaded instance ptid(7)");
    test_value(cube1.dir(11,9).dir().ra_deg(), 83.77876, "Loaded instance RA");
    test_value(cube1.dir(11,9).dir().dec_deg(), 22.03696, "Loaded instance RA");
    test_value(cube1.spi_x(11).ra_deg(), 83.77876, "Loaded instance RA_X");
    test_value(cube1.spi_x(11).dec_deg(), 22.03696, "Loaded instance DEC_X");
    test_value(cube1.spi_z(11).ra_deg(), 354.8671, "Loaded instance RA_Z");
    test_value(cube1.spi_z(11).dec_deg(), -2.686594, "Loaded instance DEC_Z");

    // Test one event bin
    const GSPIEventBin* bin1 = cube1[45];
    test_value(bin1->dir().dir().ra_deg(), 87.98045, "Bin[45] RA");
    test_value(bin1->dir().dir().dec_deg(), 17.98266, "Bin[45] DEC");
    test_value(bin1->dir().detid(), 4, "Bin[45] detector identifier");
    test_value(bin1->energy().keV(), 262.5, "Bin[45] energy");
    test_value(bin1->time().secs(), -216415756.731, "Bin[45] time");
    test_value(bin1->counts(), 4766.0, "Bin[45] counts");
    test_value(bin1->error(), 69.03622, "Bin[45] error");
    test_value(bin1->ontime(), 2928.9207878, "Bin[45] ontime");
    test_value(bin1->livetime(), 2578.60894906, "Bin[45] livetime");
    test_value(bin1->index(), 45, "Bin[45] index");
    test_value(bin1->ipt(), 0, "Bin[45] pointing index");
    test_value(bin1->idir(), 4, "Bin[45] direction index");
    test_value(bin1->iebin(), 5, "Bin[45] energy bin index");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test GSPIModelDataSpace class
 ***************************************************************************/
void TestGSPI::test_modeldataspace(void)
{
    // Test content of void model
    GSPIModelDataSpace model1;
    test_value(model1.classname(), "GSPIModelDataSpace", "Void instance classname");
    test_value(model1.type(), "DataSpace", "Void instance type");
    test_value(model1.size(), 0, "Void instance size");
    test_value(model1.instruments(), "SPI", "Void instance instruments");
    test_assert(model1.is_constant(), "Void instance is_constant");
    GXmlElement xml;
    model1.write(xml);
    model1.read(*xml.element(0));
    test_value(model1.type(), "DataSpace", "Void instance type");
    test_value(model1.size(), 0, "Void instance size");
    test_value(model1.instruments(), "SPI", "Void instance instruments");
    test_assert(model1.is_constant(), "Void instance is_constant");

    // Construct "orbit" model from OG
    GSPIObservation obs(og_dol);
    GSPIModelDataSpace model2(obs, "GEDSAT", "orbit", 0);
    test_value(model2.type(), "DataSpace", "ORBIT instance type");
    test_value(model2.size(), 1, "ORBIT instance size");
    test_value(model2.instruments(), "SPI", "ORBIT instance instruments");
    test_value(model2.name(), "GEDSAT", "ORBIT instance name");
    test_value(model2[0].name(), "GEDSAT O0044", "ORBIT instance parameter");
    xml.clear();
    model2.write(xml);

    // Construct "orbit" model from XML
    GSPIModelDataSpace model3(*xml.element(0));
    test_value(model3.type(), "DataSpace", "XML constructed instance type");
    test_value(model3.size(), 1, "XML constructed ORBIT instance size");
    test_value(model3.instruments(), "SPI", "XML constructed ORBIT instance instruments");
    test_value(model3.name(), "GEDSAT", "ORBIT instance name");
    test_value(model3[0].name(), "GEDSAT O0044", "XML constructed ORBIT instance parameter");

    // Construct "orbit,dete" model from OG
    GSPIModelDataSpace model4(obs, "GEDSAT", "orbit,dete", 0);
    test_value(model4.size(), 19, "ORBIT,DETE instance size");
    test_value(model4[1].name(), "GEDSAT D001 O0044", "ORBIT,DETE instance parameter");

    // Construct "orbit,dete,ebin" model from OG
    GSPIModelDataSpace model5(obs, "GEDSAT", "orbit,dete,ebin", 0);
    test_value(model5.size(), 190, "ORBIT,DETE,EBIN instance size");
    test_value(model5[35].name(), "GEDSAT E005 D003 O0044", "ORBIT,DETE,EBIN instance parameter");

    // Construct "point" model from OG
    GSPIModelDataSpace model6(obs, "GEDSAT", "point", 0);
    test_value(model6.size(), 88, "POINT instance size");
    test_value(model6[35].name(), "GEDSAT P000035", "POINT instance parameter");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test GSPIObservation class
 ***************************************************************************/
void TestGSPI::test_obs(void)
{
    // Test content of void observation
    GSPIObservation obs1;
    test_value(obs1.classname(), "GSPIObservation", "Void instance classname");
    test_value(obs1.instrument(), "SPI", "Void instance instrument");
    test_value(obs1.ontime(), 0.0, "Void instance ontime");
    test_value(obs1.livetime(), 0.0, "Void instance livetime");
    test_value(obs1.deadc(), 0.0, "Void instance deadc");

    // Test content of OG
    GSPIObservation obs2(og_dol);
    test_value(obs2.instrument(), "SPI", "Void instance instrument");
    test_value(obs2.ontime(), 193966.8178673, "Void instance ontime");
    test_value(obs2.livetime(), 170657.5371606, "Void instance livetime");
    test_value(obs2.deadc(), 0.87982851416, "Void instance deadc");

    // Return
    return;
}


/***********************************************************************//**
 * @brief Test GSPIResponse class
 ***************************************************************************/
void TestGSPI::test_response(void)
{
    // Set energy
    GEnergy obsEnergy;

    // Test content of void response
    GSPIResponse rsp1;
    test_value(rsp1.classname(), "GSPIResponse", "Void instance classname");
    test_assert(!rsp1.use_edisp(), "Void instance use_edisp");
    test_assert(!rsp1.use_tdisp(), "Void instance use_edisp");
    test_assert(!rsp1.is_precomputed(), "Void instance is_precomputed");
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Test GSPITools module
 ***************************************************************************/
void TestGSPI::test_tools(void)
{
    // Test methods
    test_value(gammalib::spi_ijd2time(1500.0).utc(), "2004-02-08T23:58:56",
               "spi_ijd2time");

    // Return
    return;
}


/***************************************************************************
 * @brief Main entry point for test executable
 ***************************************************************************/
int main(void)
{
    // Allocate test suit container
    GTestSuites testsuites("INTEGRAL/SPI instrument specific class testing");

    // Set CALDB environment variable
    std::string caldb = "CALDB="+spi_caldb;
    putenv((char*)caldb.c_str());

    // Initially assume that we pass all tests
    bool success = true;

    // Create test suites and append them to the container
    TestGSPI suite;
    testsuites.append(suite);

    // Run the testsuites
    success = testsuites.run();

    // Save test report
    testsuites.save("reports/GSPI.xml");

    // Return success status
    return (success ? 0 : 1);
}
