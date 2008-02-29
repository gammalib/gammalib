/***************************************************************************
 *                  test_GFits.cpp  -  test FITS classes                   *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2008 by Jurgen Knodlseder                   *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

/* __ Includes ___________________________________________________________ */
#include "test_GFits.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Globals ____________________________________________________________ */


/***************************************************************************
 *                         Verification routines                           *
 ***************************************************************************/
int fequal(double val, double ref)
{
    return (fabs(val-ref) < 1.0e-5) ? 1 : 0;
}
int dequal(double val, double ref)
{
    return (fabs(val-ref) < 1.0e-10) ? 1 : 0;
}


/***************************************************************************
 *                              Test: Open                                 *
 ***************************************************************************/
void test_open(void)
{
    std::cout << "Test GFits: Open FITS file: ";
    try {
        GFits fits;
        fits.open("test_GFits.fits");
    }
    catch (exception &e) {
        std::cout << std::endl << "TEST ERROR: Unable to open FITS file." << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ". ok." << std::endl;
}


/***************************************************************************
 *                             Test: Columns                               *
 ***************************************************************************/
void test_columns(void)
{
    GFits fits;
    try {
        // Open FITS file
        fits.open("test_GFits.fits");
    }
    catch (exception &e) {
        std::cout << std::endl << "TEST ERROR: Unable to open test file." << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }


    //
    // Test single string column
    //
    std::cout << "Test GFits: String column handling: ";
    try {
        GFitsTableStrCol str = *((GFitsTableStrCol*)fits.hdu("BinTable")->column("TSTRING"));
        if (str.string(0) != "1" || str.string(1) != "2" || str.string(2) != "") {
          cout << endl << "TEST ERROR: Bad single TSTRING values read (" 
               << str.string(0) << ", " << str.string(1) << ", " << str.string(2) << ")." << endl;
          throw;
        }
        else
          std::cout << ".";
        std::string nullstr = "empty";
        str.set_nullstr(nullstr);
        if (str.string(2) != "empty") {
          cout << endl << "TEST ERROR: Bad single TSTRING NULL value read (" 
               << str.string(2) << ")." << endl;
          throw;
        }
        else
          std::cout << ".";
    }
    catch (exception &e) {
        std::cout << std::endl << "TEST ERROR: Unable to handle single TSTRING column." << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ". ok." << std::endl;

    //
    // Test single short column
    //
    std::cout << "Test GFits: Short column handling: ";
    try {
        GFitsTableShtCol sht = *((GFitsTableShtCol*)fits.hdu("BinTable")->column("TSHORT"));
        if (sht.integer(0) !=  1 || sht.real(0) !=  1.0 || sht.string(0) !=  "1" ||
            sht.integer(1) !=  2 || sht.real(1) !=  2.0 || sht.string(1) !=  "2" ||
            sht.integer(2) != -3 || sht.real(2) != -3.0 || sht.string(2) != "-3") {
          cout << endl << "TEST ERROR: Bad single TSHORT values read." << endl;
          throw;
        }
        else
          std::cout << ".";
    }
    catch (exception &e) {
        std::cout << std::endl << "TEST ERROR: Unable to handle single TSHORT column." << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ". ok." << std::endl;

    //
    // Test single long column
    //
    std::cout << "Test GFits: Long column handling: ";
    try {
        GFitsTableLngCol lng = *((GFitsTableLngCol*)fits.hdu("BinTable")->column("TLONG"));
        if (lng.integer(0) !=  10 || lng.real(0) !=  10.0 || lng.string(0) !=  "10" ||
            lng.integer(1) != -20 || lng.real(1) != -20.0 || lng.string(1) != "-20" ||
            lng.integer(2) !=  30 || lng.real(2) !=  30.0 || lng.string(2) !=  "30") {
          cout << endl << "TEST ERROR: Bad single TLONG values read." << endl;
          throw;
        }
        else
          std::cout << ".";
    }
    catch (exception &e) {
        std::cout << std::endl << "TEST ERROR: Unable to handle single TLONG column." << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ". ok." << std::endl;

    //
    // Test single floating point column
    //
    std::cout << "Test GFits: Floating point column handling: ";
    try {
        GFitsTableFltCol flt = *((GFitsTableFltCol*)fits.hdu("BinTable")->column("TFLOAT"));
        if (!fequal(flt.real(0),1.0)) {
          cout << endl << "TEST ERROR: Bad single TFLOAT values read (" 
               << flt.real(0) << " != 1.0)." << endl;
          throw;
        }
        else
          std::cout << ".";
        if (!fequal(flt.real(1),-2.1)) {
          cout << endl << "TEST ERROR: Bad single TFLOAT values read (" 
               << flt.real(1) << " != -2.1)." << endl;
          throw;
        }
        else
          std::cout << ".";
        if (flt.integer(0) != 1) {
          cout << endl << "TEST ERROR: Bad single TFLOAT values read (" 
               << flt.integer(0) << " != 1)." << endl;
          throw;
        }
        else
          std::cout << ".";
        if (flt.integer(1) != -2) {
          cout << endl << "TEST ERROR: Bad single TFLOAT values read (" 
               << flt.integer(1) << " != -2)." << endl;
          throw;
        }
        else
          std::cout << ".";
        if (flt.string(0) != "1.000000e+00") {
          cout << endl << "TEST ERROR: Bad single TFLOAT values read (" 
               << flt.string(0) << " != 1.000000e+00)." << endl;
          throw;
        }
        else
          std::cout << ".";
        if (flt.string(1) != "-2.100000e+00") {
          cout << endl << "TEST ERROR: Bad single TFLOAT values read (" 
               << flt.string(1) << " != -2.100000e+00)." << endl;
          throw;
        }
        else
          std::cout << ".";
        float nullval = -99.9;
        flt.set_nullval(&nullval);
        if (!fequal(flt.real(2),-99.9)) {
          cout << endl << "TEST ERROR: Bad single TFLOAT NULL value read (" 
               << flt.real(2) << ")." << endl;
          throw;
        }
        else
          std::cout << ".";
    }
    catch (exception &e) {
        std::cout << std::endl << "TEST ERROR: Unable to handle single TFLOAT column." << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    try {
        GFitsTableFltCol flt = *((GFitsTableFltCol*)fits.hdu("BinTable2")->column("TFLOAT3"));
        if (!fequal(flt.real(0,0),0.0)) {
          cout << endl << "TEST ERROR: Bad TFLOAT (3) values read (" 
               << flt.real(0,0) << " != 0.0)." << endl;
          throw;
        }
    }
    catch (exception &e) {
        std::cout << std::endl << "TEST ERROR: Unable to handle single TFLOAT column." << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ". ok." << std::endl;

    //
    // Test single double column
    //
    std::cout << "Test GFits: Double precision column handling: ";
    try {
        GFitsTableDblCol dbl = *((GFitsTableDblCol*)fits.hdu("BinTable")->column("TDOUBLE"));
        if (!dequal(dbl.real(0),11.1)) {
          cout << endl << "TEST ERROR: Bad single TDOUBLE values read (" 
               << dbl.real(0) << " != 11.1)." << endl;
          throw;
        }
        else
          std::cout << ".";
        if (!dequal(dbl.real(2),22.2)) {
          cout << endl << "TEST ERROR: Bad single TDOUBLE values read (" 
               << dbl.real(2) << " != 22.2)." << endl;
          throw;
        }
        else
          std::cout << ".";
        if (dbl.integer(0) != 11) {
          cout << endl << "TEST ERROR: Bad single TDOUBLE values read (" 
               << dbl.integer(0) << " != 11)." << endl;
          throw;
        }
        else
          std::cout << ".";
        if (dbl.integer(2) != 22) {
          cout << endl << "TEST ERROR: Bad single TDOUBLE values read (" 
               << dbl.integer(2) << " != 22)." << endl;
          throw;
        }
        else
          std::cout << ".";
        if (dbl.string(0) != "1.110000e+01") {
          cout << endl << "TEST ERROR: Bad single TDOUBLE values read (" 
               << dbl.string(0) << " != 1.110000e+01)." << endl;
          throw;
        }
        else
          std::cout << ".";
        if (dbl.string(2) != "2.220000e+01") {
          cout << endl << "TEST ERROR: Bad single TDOUBLE values read (" 
               << dbl.string(2) << " != 2.220000e+01)." << endl;
          throw;
        }
        else
          std::cout << ".";
        double nullval = -9999.9;
        dbl.set_nullval(&nullval);
        if (!dequal(dbl.real(1),-9999.9)) {
          cout << endl << "TEST ERROR: Bad single TDOUBLE NULL value read (" 
               << dbl.real(1) << ")." << endl;
          throw;
        }
        else
          std::cout << ".";
    }
    catch (exception &e) {
        std::cout << std::endl << "TEST ERROR: Unable to handle single TDOUBLE column." << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ". ok." << std::endl;

}


/***************************************************************************
 *                            Main test function                           *
 ***************************************************************************/
int main(void)
{
    // Dump header
    std::cout << std::endl;
    std::cout << "***********************" << std::endl;
    std::cout << "* GFits class testing *" << std::endl;
    std::cout << "***********************" << std::endl;

    // Execute the tests
    test_open();
    test_columns();

    // Return
    return 0;
}
