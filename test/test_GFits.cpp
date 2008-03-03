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
#include <stdlib.h>
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
    cout << "Test GFits: Open FITS file: ";
    try {
        GFits fits;
        fits.open("test_GFits.fits");
        cout << endl << fits;
    }
    catch (exception &e) {
        cout << endl << "TEST ERROR: Unable to open FITS file." << endl;
        cout << e.what() << endl;
        throw;
    }
    cout << ". ok." << endl;
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
        cout << endl << "TEST ERROR: Unable to open test file." << endl;
        cout << e.what() << endl;
        throw;
    }

    //
    // Test single string column
    //
    cout << "Test GFits: String column handling: ";
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
        cout << endl << "TEST ERROR: Unable to handle single TSTRING column." << endl;
        cout << e.what() << endl;
        throw;
    }
    cout << ". ok." << endl;

    //
    // Test single short column
    //
    cout << "Test GFits: Short column handling: ";
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
        cout << endl << "TEST ERROR: Unable to handle single TSHORT column." << endl;
        cout << e.what() << endl;
        throw;
    }
    cout << ". ok." << endl;

    //
    // Test single long column
    //
    cout << "Test GFits: Long column handling: ";
    try {
        GFitsTableLngCol lng = *((GFitsTableLngCol*)fits.hdu("BinTable")->column("TLONG"));
        if (lng.integer(0) !=  10 || lng.real(0) !=  10.0 || lng.string(0) !=  "10" ||
            lng.integer(1) != -20 || lng.real(1) != -20.0 || lng.string(1) != "-20" ||
            lng.integer(2) !=  30 || lng.real(2) !=  30.0 || lng.string(2) !=  "30") {
          cout << endl << "TEST ERROR: Bad single TLONG values read." << endl;
          throw;
        }
        else
          cout << ".";
    }
    catch (exception &e) {
        cout << endl << "TEST ERROR: Unable to handle single TLONG column." << endl;
        cout << e.what() << endl;
        throw;
    }
    cout << ". ok." << endl;

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
        cout << endl << "TEST ERROR: Unable to handle single TFLOAT column." << endl;
        cout << e.what() << endl;
        throw;
    }
    cout << ". ok." << endl;

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
        cout << endl << "TEST ERROR: Unable to handle single TDOUBLE column." << endl;
        cout << e.what() << endl;
        throw;
    }
    cout << ". ok." << endl;

}


/***************************************************************************
 *                            Test: Creation                               *
 ***************************************************************************/
void test_create(void)
{
    // Dump header
    cout << "Test GFits: Create FITS file: ";

    // Remove FITS file
    system("rm -rf test.fits");

    // Create FITS file with empty double precision image
    try {
        GFits fits;
        fits.open("test.fits");
        GFitsDblImage image;
        GFitsHDU hdu(image);
        fits.append(&hdu);
        //cout << fits;
        fits.save();
    }
    catch (exception &e) {
        cout << endl << "TEST ERROR: Unable to create FITS file." << endl;
        cout << e.what() << endl;
        throw;
    }
    cout << ".";

    // Attach double precision image
    double sum = 0.0;
    try {
        GFits fits;
        fits.open("test.fits");
        int naxis       = 2;
        int nx          = 10;
        int ny          = 20;
        int naxes[]     = {nx,ny};
        GFitsDblImage image(naxis, naxes);
        for (int ix = 0; ix < nx; ++ix) {
            for (int iy = 0; iy < ny; ++iy) {
                image(ix,iy) = 0.01*ix + 0.01*iy;
                sum += image(ix,iy);
            }
        }
        GFitsHDU hdu(image);
        fits.append(&hdu);
        //cout << fits;
        fits.save();
    }
    catch (exception &e) {
        cout << endl << "TEST ERROR: Unable to create FITS file." << endl;
        cout << e.what() << endl;
        throw;
    }
    cout << ".";

    // Re-open double precision image
    try {
        GFits fits;
        fits.open("test.fits");
        GFitsHDU*      hdu   = fits.hdu(2);
        GFitsDblImage* image = (GFitsDblImage*)hdu->data();
        int nx = image->naxes(0);
        int ny = image->naxes(1);
        double total = 0.0;
        for (int ix = 0; ix < nx; ++ix) {
            for (int iy = 0; iy < ny; ++iy) {
                total += (*image)(ix,iy);
            }
        }
        if (!dequal(total, sum)) {
            cout << endl << "TEST ERROR: Bad values in loaded image." << endl;
            throw;
        }
    }
    catch (exception &e) {
        cout << endl << "TEST ERROR: Unable to re-open FITS file." << endl;
        cout << e.what() << endl;
        throw;
    }
    cout << ".";

    // Attach binary table
    try {
        GFits fits;
        fits.open("test.fits");
        int nrows = 10;
        int ncols = 5;
        GFitsBinTable table(nrows,ncols);
        GFitsHDU hdu(table);
        fits.append(&hdu);
        cout << fits;
        fits.save();
    }
    catch (exception &e) {
        cout << endl << "TEST ERROR: Unable to attach binary table." << endl;
        cout << e.what() << endl;
        throw;
    }
    cout << ".";

    // Attach ASCII table
    try {
        GFits fits;
        fits.open("test.fits");
        int nrows = 10;
        int ncols = 5;
        GFitsAsciiTable table(nrows,ncols);
        GFitsHDU hdu(table);
        fits.append(&hdu);
        cout << fits;
        fits.save();
    }
    catch (exception &e) {
        cout << endl << "TEST ERROR: Unable to attach ASCII table." << endl;
        cout << e.what() << endl;
        throw;
    }
    cout << ".";
        

    cout << ". ok." << endl;
}


/***************************************************************************
 *                        Test: GFitsDblImage test                         *
 ***************************************************************************/
void test_image_double(void)
{
    // Dump header
    cout << "Test GFitsDblImage: ";
    
    // Test 1D Image
    try {
        int naxis       = 1;
        int nx          = 10;
        int naxes[]     = {nx};
        GFitsDblImage image(naxis, naxes);
        cout << ".";
        
        // Fill and read image
        try {
            for (int ix = 0; ix < nx; ++ix)
                image(ix) = double(ix);
            double sum = 0.0;
            for (int ix = 0; ix < nx; ++ix)
                sum += image(ix);
            if (!dequal(sum, 0.5*double(nx-1)*double(nx))) {
                cout << endl << "TEST ERROR: Bad fill/read of 1D image." << endl;
                throw;
            }
        }
        catch (exception &e) {
            cout << endl << "TEST ERROR: Unable to fill/read 1D image." << endl;
            cout << e.what() << endl;
            throw;
        }
        cout << ".";

        // Bad 2D operator
        try {
            double result = image(0,1);
        }
        catch (GException::fits_wrong_image_operator &e) {
            cout << ".";
        }
        catch (exception &e) {
            cout << endl << "TEST ERROR: Wrong operactor access error." << endl;
            cout << e.what() << endl;
            throw;
        }

        // Bad 3D operator
        try {
            double result = image(0,1,2);
        }
        catch (GException::fits_wrong_image_operator &e) {
            cout << ".";
        }
        catch (exception &e) {
            cout << endl << "TEST ERROR: Wrong operactor access error." << endl;
            cout << e.what() << endl;
            throw;
        }

        // Bad 4D operator
        try {
            double result = image(0,1,2,3);
        }
        catch (GException::fits_wrong_image_operator &e) {
            cout << ".";
        }
        catch (exception &e) {
            cout << endl << "TEST ERROR: Wrong operactor access error." << endl;
            cout << e.what() << endl;
            throw;
        }

    }
    catch (exception &e) {
        cout << endl << "TEST ERROR: Unable to test GFitsDblImage classes." << endl;
        cout << e.what() << endl;
        throw;
    }
    cout << ". ok." << endl;
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
    test_image_double();
    test_create();
    //test_open();
    //test_columns();

    // Return
    return 0;
}
