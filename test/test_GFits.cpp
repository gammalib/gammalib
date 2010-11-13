/***************************************************************************
 *                  test_GFits.cpp  -  test FITS classes                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2010 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

/* __ Includes ___________________________________________________________ */
#include <iostream>                           // cout, cerr
#include <stdexcept>                          // std::exception
#include <stdlib.h>
#include "GammaLib.hpp"

/* __ Namespaces _________________________________________________________ */
using namespace std;

/* __ Globals ____________________________________________________________ */


/***************************************************************************
 *                         Verification routines                           *
 ***************************************************************************/
int equal(double val, double ref, double eps)
{
    return (fabs(val-ref) < eps) ? 1 : 0;
}
int fequal(double val, double ref)
{
    return (fabs(val-ref) < 1.0e-5) ? 1 : 0;
}
int dequal(double val, double ref)
{
    return (fabs(val-ref) < 1.0e-10) ? 1 : 0;
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
        GFitsImageDbl image;
        GFitsHDU hdu(image);
        fits.append(hdu);
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
        GFitsImageDbl image(naxis, naxes);
        for (int ix = 0; ix < nx; ++ix) {
            for (int iy = 0; iy < ny; ++iy) {
                image(ix,iy) = 0.01*ix + 0.01*iy;
                sum += image(ix,iy);
            }
        }
        GFitsHDU hdu(image);
        fits.append(hdu);
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
        GFitsHDU*      hdu   = fits.hdu(1);
        GFitsImageDbl* image = (GFitsImageDbl*)hdu->data();
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
        // Re-open FITS file
        GFits fits;
        fits.open("test.fits");

        // Create binary Table with 10 rows
        int nrows = 10;
        GFitsBinTable    table  = GFitsBinTable(nrows);
        GFitsTableDblCol first  = GFitsTableDblCol("First", nrows);
        GFitsTableDblCol second = GFitsTableDblCol("Second", nrows);
        first(0)  =  1.0;
        first(1)  =  2.0;
        first(2)  =  3.0;
        second(0) = -99.0;
        cout << endl << first << endl;
        cout << second << endl;
        table.append_column(first);
        table.append_column(second);
//        cout << table;

        // Create HDU and append to FILE file
        GFitsHDU hdu(table);
        fits.append(hdu);

        // Save FITS file
        fits.save();
    }
    catch (exception &e) {
        cout << endl << "TEST ERROR: Unable to attach binary table." << endl;
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
        GFitsImageDbl image(naxis, naxes);
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
 *                  Test: double precision FITS binary table               *
 ***************************************************************************/
void test_bintable_double(void)
{
    // Set filename
    std::string filename = "test_bintable_double.fits";

    // Dump header
    cout << "Test GFitsTableDblCol: ";
    
    // Remove FITS file
    std::string cmd = "rm -rf "+ filename;
    system(cmd.c_str());

    // Allocate reference sums
    double      sum_dbl       = 0.0;
    double      sum_dbl10     = 0.0;
    int         sum_dbl_int   = 0;
    int         sum_dbl10_int = 0;
    std::string sum_dbl_str;
    std::string sum_dbl10_str;

    // Build tables
    try {
        // Open FITS file
        GFits fits;
        fits.open(filename);

        // Set number of rows
        int nrows = 10;
        int nvec  = 10;

        // Initial control sums
        float       tot_flt = 0.0;
        double      tot_dbl = 0.0;
        int         tot_int = 0;
        long        tot_lng = 0;
        short       tot_sht = 0;
        std::string tot_str;

        //
        // ===== D O U B L E =====
        //

        // Set double precision table
        GFitsTableDblCol col_dbl = GFitsTableDblCol("DOUBLE", nrows);
        for (int i = 0; i < nrows; ++i) {
            double val_dbl = cos(double(i));
            col_dbl(i)   = val_dbl;
            sum_dbl     += val_dbl;
            sum_dbl_int += int(val_dbl);
            ostringstream s_value;
            s_value << scientific << val_dbl;
            sum_dbl_str += ":"+s_value.str();
        }
        cout << ".";

        // Check double precision table (operator access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i)
            tot_dbl += col_dbl(i);
        if (!dequal(tot_dbl, sum_dbl)) {
            cout << endl << "TEST ERROR: GFitsTableDblCol::operator()";
            cout << endl << "  Reference sum: " << sum_dbl;
            cout << endl << "  Derived sum:   " << tot_dbl << endl;
            throw;
        }
        cout << ".";

        // Check double precision table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i)
            tot_dbl += col_dbl.real(i);
        if (!dequal(tot_dbl, sum_dbl)) {
            cout << endl << "TEST ERROR: GFitsTableDblCol::real()";
            cout << endl << "  Reference sum: " << sum_dbl;
            cout << endl << "  Derived sum:   " << tot_dbl << endl;
            throw;
        }
        cout << ".";

        // Check double precision table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i)
            tot_int += col_dbl.integer(i);
        if (tot_int != sum_dbl_int) {
            cout << endl << "TEST ERROR: GFitsTableDblCol::integer()";
            cout << endl << "  Reference sum: " << sum_dbl_int;
            cout << endl << "  Derived sum:   " << tot_int << endl;
            throw;
        }
        cout << ".";

        // Check double precision table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i)
            tot_str += ":"+col_dbl.string(i);
        if (tot_str != sum_dbl_str) {
            cout << endl << "TEST ERROR: GFitsTableDblCol::string()";
            cout << endl << "  Reference sum: " << sum_dbl_str;
            cout << endl << "  Derived sum:   " << tot_str << endl;
            throw;
        }
        cout << ".";

        //
        // ===== D O U B L E 1 0 =====
        //

        // Set double precision table
        GFitsTableDblCol col_dbl10 = GFitsTableDblCol("DOUBLE10", nrows, nvec);
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j) {
                double val_dbl = cos(double(i))*cos(0.33*double(j));
                col_dbl10(i,j) = val_dbl;
                sum_dbl10     += val_dbl;
                sum_dbl10_int += int(val_dbl);
                ostringstream s_value;
                s_value << scientific << val_dbl;
                sum_dbl10_str += ":"+s_value.str();
            }
        }
        cout << ".";

        // Check double precision table (operator access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_dbl += col_dbl10(i,j);
        }
        if (!dequal(tot_dbl, sum_dbl10)) {
            cout << endl << "TEST ERROR: GFitsTableDblCol::operator() - 10";
            cout << endl << "  Reference sum: " << sum_dbl10;
            cout << endl << "  Derived sum:   " << tot_dbl << endl;
            throw;
        }
        cout << ".";

        // Check double precision table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_dbl += col_dbl10.real(i,j);
        }
        if (!dequal(tot_dbl, sum_dbl10)) {
            cout << endl << "TEST ERROR: GFitsTableDblCol::real() - 10";
            cout << endl << "  Reference sum: " << sum_dbl10;
            cout << endl << "  Derived sum:   " << tot_dbl << endl;
            throw;
        }
        cout << ".";

        // Check double precision table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_int += col_dbl10.integer(i,j);
        }
        if (tot_int != sum_dbl10_int) {
            cout << endl << "TEST ERROR: GFitsTableDblCol::integer() - 10";
            cout << endl << "  Reference sum: " << sum_dbl10_int;
            cout << endl << "  Derived sum:   " << tot_int << endl;
            throw;
        }
        cout << ".";

        // Check double precision table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_str += ":"+col_dbl10.string(i,j);
        }
        if (tot_str != sum_dbl10_str) {
            cout << endl << "TEST ERROR: GFitsTableDblCol::string() - 10";
            cout << endl << "  Reference sum: " << sum_dbl10_str;
            cout << endl << "  Derived sum:   " << tot_str << endl;
            throw;
        }
        cout << ".";

        //
        // ===== W R I T E   T A B L E =====
        //

        // Build binary table
        GFitsBinTable table = GFitsBinTable(nrows);
        table.append_column(col_dbl);
        table.insert_column(1, col_dbl10);

        // Create HDU and append to FILE file
        GFitsHDU hdu(table);
        fits.append(hdu);

        // Save FITS file
        fits.save();
        cout << ".";
    }
    catch (exception &e) {
        cout << endl << "TEST ERROR: Unable to build tables." << endl;
        cout << e.what() << endl;
        throw;
    }
    cout << ".";


    //================================================================
    //================================================================
    //================================================================
    // Read tables
    //================================================================
    //================================================================
    //================================================================
    try {
        // Open FITS file
        GFits fits;
        fits.open(filename);

        // Set number of rows
        int nrows = 10;
        int nvec  = 10;

        // Initial control sums
        float       tot_flt = 0.0;
        double      tot_dbl = 0.0;
        int         tot_int = 0;
        std::string tot_str;

        //
        // ===== D O U B L E =====
        //

        // Get column
        GFitsTableDblCol* col_dbl =
                    (GFitsTableDblCol*)fits.hdu(1)->column("DOUBLE");

        // Check double precision table (operator access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i)
            tot_dbl += (*col_dbl)(i);
        if (!dequal(tot_dbl, sum_dbl)) {
            cout << endl << "TEST ERROR: GFitsTableDblCol::operator()";
            cout << endl << "  Reference sum: " << sum_dbl;
            cout << endl << "  Derived sum:   " << tot_dbl << endl;
            throw;
        }
        cout << ".";

        // Check double precision table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i)
            tot_dbl += (*col_dbl).real(i);
        if (!dequal(tot_dbl, sum_dbl)) {
            cout << endl << "TEST ERROR: GFitsTableDblCol::real()";
            cout << endl << "  Reference sum: " << sum_dbl;
            cout << endl << "  Derived sum:   " << tot_dbl << endl;
            throw;
        }
        cout << ".";

        // Check double precision table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i)
            tot_int += (*col_dbl).integer(i);
        if (tot_int != sum_dbl_int) {
            cout << endl << "TEST ERROR: GFitsTableDblCol::integer()";
            cout << endl << "  Reference sum: " << sum_dbl_int;
            cout << endl << "  Derived sum:   " << tot_int << endl;
            throw;
        }
        cout << ".";

        // Check double precision table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i)
            tot_str += ":"+ (*col_dbl).string(i);
        if (tot_str != sum_dbl_str) {
            cout << endl << "TEST ERROR: GFitsTableDblCol::string()";
            cout << endl << "  Reference sum: " << sum_dbl_str;
            cout << endl << "  Derived sum:   " << tot_str << endl;
            throw;
        }
        cout << ".";

        //
        // ===== D O U B L E 1 0 =====
        //

        // Get column
        GFitsTableDblCol* col_dbl10 =
                    (GFitsTableDblCol*)fits.hdu(1)->column("DOUBLE10");

        // Check double precision table (operator access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_dbl += (*col_dbl10)(i,j);
        }
        if (!dequal(tot_dbl, sum_dbl10)) {
            cout << endl << "TEST ERROR: GFitsTableDblCol::operator() - 10";
            cout << endl << "  Reference sum: " << sum_dbl10;
            cout << endl << "  Derived sum:   " << tot_dbl << endl;
            throw;
        }
        cout << ".";

        // Check double precision table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_dbl += col_dbl10->real(i,j);
        }
        if (!dequal(tot_dbl, sum_dbl10)) {
            cout << endl << "TEST ERROR: GFitsTableDblCol::real() - 10";
            cout << endl << "  Reference sum: " << sum_dbl10;
            cout << endl << "  Derived sum:   " << tot_dbl << endl;
            throw;
        }
        cout << ".";

        // Check double precision table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_int += col_dbl10->integer(i,j);
        }
        if (tot_int != sum_dbl10_int) {
            cout << endl << "TEST ERROR: GFitsTableDblCol::integer() - 10";
            cout << endl << "  Reference sum: " << sum_dbl10_int;
            cout << endl << "  Derived sum:   " << tot_int << endl;
            throw;
        }
        cout << ".";

        // Check double precision table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_str += ":"+col_dbl10->string(i,j);
        }
        if (tot_str != sum_dbl10_str) {
            cout << endl << "TEST ERROR: GFitsTableDblCol::string() - 10";
            cout << endl << "  Reference sum: " << sum_dbl10_str;
            cout << endl << "  Derived sum:   " << tot_str << endl;
            throw;
        }
        cout << ".";
    }
    catch (exception &e) {
        cout << endl << "TEST ERROR: Unable to build tables." << endl;
        cout << e.what() << endl;
        throw;
    }
    cout << ". ok." << endl;

}


/***************************************************************************
 *                  Test: single precision FITS binary table               *
 ***************************************************************************/
void test_bintable_float(void)
{
    // Set filename
    std::string filename = "test_bintable_float.fits";

    // Dump header
    cout << "Test GFitsTableFltCol: ";
    
    // Remove FITS file
    std::string cmd = "rm -rf "+ filename;
    system(cmd.c_str());

    // Allocate reference sums
    float       sum_flt       = 0.0;
    float       sum_flt10     = 0.0;
    int         sum_flt_int   = 0;
    int         sum_flt10_int = 0;
    std::string sum_flt_str;
    std::string sum_flt10_str;

    // Build tables
    try {
        // Open FITS file
        GFits fits;
        fits.open(filename);

        // Set number of rows
        int nrows = 10;
        int nvec  = 10;

        // Initial control sums
        float       tot_flt = 0.0;
        double      tot_dbl = 0.0;
        int         tot_int = 0;
        long        tot_lng = 0;
        short       tot_sht = 0;
        std::string tot_str;

        //
        // ===== F L O A T =====
        //

        // Set single precision table
        GFitsTableFltCol col_flt = GFitsTableFltCol("FLOAT", nrows);
        for (int i = 0; i < nrows; ++i) {
            float val_flt = cos(0.1*float(i));
            col_flt(i)   = val_flt;
            sum_flt     += val_flt;
            sum_flt_int += int(val_flt);
            ostringstream s_value;
            s_value << scientific << val_flt;
            sum_flt_str += ":"+s_value.str();
        }
        cout << ".";

        // Check single precision table (operator access)
        tot_flt = 0.0;
        for (int i = 0; i < nrows; ++i)
            tot_flt += col_flt(i);
        if (!fequal(tot_flt, sum_flt)) {
            cout << endl << "TEST ERROR: GFitsTableFltCol::operator()";
            cout << endl << "  Reference sum: " << sum_flt;
            cout << endl << "  Derived sum:   " << tot_flt << endl;
            throw;
        }
        cout << ".";

        // Check single precision table (real access)
        tot_flt = 0.0;
        for (int i = 0; i < nrows; ++i)
            tot_flt += col_flt.real(i);
        if (!fequal(tot_flt, sum_flt)) {
            cout << endl << "TEST ERROR: GFitsTableFltCol::real()";
            cout << endl << "  Reference sum: " << sum_flt;
            cout << endl << "  Derived sum:   " << tot_flt << endl;
            throw;
        }
        cout << ".";

        // Check single precision table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i)
            tot_int += col_flt.integer(i);
        if (tot_int != sum_flt_int) {
            cout << endl << "TEST ERROR: GFitsTableFltCol::integer()";
            cout << endl << "  Reference sum: " << sum_flt_int;
            cout << endl << "  Derived sum:   " << tot_int << endl;
            throw;
        }
        cout << ".";

        // Check single precision table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i)
            tot_str += ":"+col_flt.string(i);
        if (tot_str != sum_flt_str) {
            cout << endl << "TEST ERROR: GFitsTableFltCol::string()";
            cout << endl << "  Reference sum: " << sum_flt_str;
            cout << endl << "  Derived sum:   " << tot_str << endl;
            throw;
        }
        cout << ".";

        //
        // ===== F L O A T  1 0 =====
        //

        // Set single precision table
        GFitsTableFltCol col_flt10 = GFitsTableFltCol("FLOAT10", nrows, nvec);
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j) {
                float val_flt  = cos(0.1*float(i))*cos(0.33*float(j));
                col_flt10(i,j) = val_flt;
                sum_flt10     += val_flt;
                sum_flt10_int += int(val_flt);
                ostringstream s_value;
                s_value << scientific << val_flt;
                sum_flt10_str += ":"+s_value.str();
            }
        }
        cout << ".";

        // Check single precision table (operator access)
        tot_flt = 0.0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_flt += col_flt10(i,j);
        }
        if (!fequal(tot_flt, sum_flt10)) {
            cout << endl << "TEST ERROR: GFitsTableFltCol::operator() - 10";
            cout << endl << "  Reference sum: " << sum_flt10;
            cout << endl << "  Derived sum:   " << tot_flt << endl;
            throw;
        }
        cout << ".";

        // Check single precision table (real access)
        tot_flt = 0.0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_flt += col_flt10.real(i,j);
        }
        if (!fequal(tot_flt, sum_flt10)) {
            cout << endl << "TEST ERROR: GFitsTableFltCol::real() - 10";
            cout << endl << "  Reference sum: " << sum_flt10;
            cout << endl << "  Derived sum:   " << tot_flt << endl;
            throw;
        }
        cout << ".";

        // Check double precision table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_int += col_flt10.integer(i,j);
        }
        if (tot_int != sum_flt10_int) {
            cout << endl << "TEST ERROR: GFitsTableFltCol::integer() - 10";
            cout << endl << "  Reference sum: " << sum_flt10_int;
            cout << endl << "  Derived sum:   " << tot_int << endl;
            throw;
        }
        cout << ".";

        // Check single precision table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_str += ":"+col_flt10.string(i,j);
        }
        if (tot_str != sum_flt10_str) {
            cout << endl << "TEST ERROR: GFitsTableFltCol::string() - 10";
            cout << endl << "  Reference sum: " << sum_flt10_str;
            cout << endl << "  Derived sum:   " << tot_str << endl;
            throw;
        }
        cout << ".";

        //
        // ===== W R I T E   T A B L E =====
        //

        // Build binary table
        GFitsBinTable table = GFitsBinTable(nrows);
        table.append_column(col_flt10);
        table.insert_column(0, col_flt);

        // Create HDU and append to FILE file
        GFitsHDU hdu(table);
        fits.append(hdu);

        // Save FITS file
        fits.save();
        cout << ".";
    }
    catch (exception &e) {
        cout << endl << "TEST ERROR: Unable to test floating point column in"
                        " binary tables." << endl;
        cout << e.what() << endl;
        throw;
    }
    cout << ".";


    //================================================================
    //================================================================
    //================================================================
    // Read tables
    //================================================================
    //================================================================
    //================================================================
    try {
        // Open FITS file
        GFits fits;
        fits.open(filename);

        // Set number of rows
        int nrows = 10;
        int nvec  = 10;

        // Initial control sums
        float       tot_flt = 0.0;
        double      tot_dbl = 0.0;
        int         tot_int = 0;
        std::string tot_str;

        //
        // ===== F L O A T =====
        //

        // Get column
        GFitsTableFltCol* col_flt =
                    (GFitsTableFltCol*)fits.hdu(1)->column("FLOAT");

        // Check single precision table (operator access)
        tot_flt = 0.0;
        for (int i = 0; i < nrows; ++i)
            tot_flt += (*col_flt)(i);
        if (!fequal(tot_flt, sum_flt)) {
            cout << endl << "TEST ERROR: GFitsTableFltCol::operator()";
            cout << endl << "  Reference sum: " << sum_flt;
            cout << endl << "  Derived sum:   " << tot_flt << endl;
            throw;
        }
        cout << ".";

        // Check single precision table (real access)
        tot_flt = 0.0;
        for (int i = 0; i < nrows; ++i)
            tot_flt += col_flt->real(i);
        if (!fequal(tot_flt, sum_flt)) {
            cout << endl << "TEST ERROR: GFitsTableFltCol::real()";
            cout << endl << "  Reference sum: " << sum_flt;
            cout << endl << "  Derived sum:   " << tot_flt << endl;
            throw;
        }
        cout << ".";

        // Check single precision table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i)
            tot_int += col_flt->integer(i);
        if (tot_int != sum_flt_int) {
            cout << endl << "TEST ERROR: GFitsTableFltCol::integer()";
            cout << endl << "  Reference sum: " << sum_flt_int;
            cout << endl << "  Derived sum:   " << tot_int << endl;
            throw;
        }
        cout << ".";

        // Check single precision table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i)
            tot_str += ":"+col_flt->string(i);
        if (tot_str != sum_flt_str) {
            cout << endl << "TEST ERROR: GFitsTableFltCol::string()";
            cout << endl << "  Reference sum: " << sum_flt_str;
            cout << endl << "  Derived sum:   " << tot_str << endl;
            throw;
        }
        cout << ".";

        //
        // ===== F L O A T  1 0 =====
        //

        // Get column
        GFitsTableFltCol* col_flt10 =
                    (GFitsTableFltCol*)fits.hdu(1)->column("FLOAT10");

        // Check single precision table (operator access)
        tot_flt = 0.0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_flt += (*col_flt10)(i,j);
        }
        if (!fequal(tot_flt, sum_flt10)) {
            cout << endl << "TEST ERROR: GFitsTableFltCol::operator() - 10";
            cout << endl << "  Reference sum: " << sum_flt10;
            cout << endl << "  Derived sum:   " << tot_flt << endl;
            throw;
        }
        cout << ".";

        // Check single precision table (real access)
        tot_flt = 0.0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_flt += col_flt10->real(i,j);
        }
        if (!fequal(tot_flt, sum_flt10)) {
            cout << endl << "TEST ERROR: GFitsTableFltCol::real() - 10";
            cout << endl << "  Reference sum: " << sum_flt10;
            cout << endl << "  Derived sum:   " << tot_flt << endl;
            throw;
        }
        cout << ".";

        // Check double precision table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_int += col_flt10->integer(i,j);
        }
        if (tot_int != sum_flt10_int) {
            cout << endl << "TEST ERROR: GFitsTableFltCol::integer() - 10";
            cout << endl << "  Reference sum: " << sum_flt10_int;
            cout << endl << "  Derived sum:   " << tot_int << endl;
            throw;
        }
        cout << ".";

        // Check single precision table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_str += ":"+col_flt10->string(i,j);
        }
        if (tot_str != sum_flt10_str) {
            cout << endl << "TEST ERROR: GFitsTableFltCol::string() - 10";
            cout << endl << "  Reference sum: " << sum_flt10_str;
            cout << endl << "  Derived sum:   " << tot_str << endl;
            throw;
        }
        cout << ".";

    }
    catch (exception &e) {
        cout << endl << "TEST ERROR: Unable to build tables." << endl;
        cout << e.what() << endl;
        throw;
    }
    cout << ". ok." << endl;

}


/***************************************************************************
 *                       Test: FITS binary table test                      *
 ***************************************************************************/
void test_bintable_short(void)
{
    // Set filename
    std::string filename = "test_bintable_short.fits";

    // Dump header
    cout << "Test GFitsTableShtCol: ";
    
    // Remove FITS file
    std::string cmd = "rm -rf "+ filename;
    system(cmd.c_str());

    // Allocate reference sums
    short       sum_sht       = 0;
    short       sum_sht10     = 0;
    short       sum_sht_int   = 0;
    short       sum_sht10_int = 0;
    std::string sum_sht_str;
    std::string sum_sht10_str;

    // Build tables
    try {
        // Open FITS file
        GFits fits;
        fits.open(filename);

        // Set number of rows
        int nrows = 10;
        int nvec  = 10;

        // Initial control sums
        float       tot_flt = 0.0;
        double      tot_dbl = 0.0;
        int         tot_int = 0;
        long        tot_lng = 0;
        short       tot_sht = 0;
        std::string tot_str;

        //
        // ===== S H O R T =====
        //

        // Set short table
        GFitsTableShtCol col_sht = GFitsTableShtCol("SHORT", nrows);
        for (int i = 0; i < nrows; ++i) {
            short val_sht = short(1000.0 * cos(0.1*float(i)));
            col_sht(i)    = val_sht;
            sum_sht      += val_sht;
            sum_sht_int  += int(val_sht);
            ostringstream s_value;
            s_value << val_sht;
            sum_sht_str += ":"+s_value.str();
        }
        cout << ".";

        // Check short table (operator access)
        tot_sht = 0;
        for (int i = 0; i < nrows; ++i)
            tot_sht += col_sht(i);
        if (tot_sht != sum_sht) {
            cout << endl << "TEST ERROR: GFitsTableShtCol::operator()";
            cout << endl << "  Reference sum: " << sum_sht;
            cout << endl << "  Derived sum:   " << tot_sht << endl;
            throw;
        }
        cout << ".";

        // Check short table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i)
            tot_dbl += col_sht.real(i);
        if (!dequal(tot_dbl, double(sum_sht))) {
            cout << endl << "TEST ERROR: GFitsTableShtCol::real()";
            cout << endl << "  Reference sum: " << sum_sht;
            cout << endl << "  Derived sum:   " << tot_dbl << endl;
            throw;
        }
        cout << ".";

        // Check short table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i)
            tot_int += col_sht.integer(i);
        if (tot_int != sum_sht_int) {
            cout << endl << "TEST ERROR: GFitsTableShtCol::integer()";
            cout << endl << "  Reference sum: " << sum_sht_int;
            cout << endl << "  Derived sum:   " << tot_int << endl;
            throw;
        }
        cout << ".";

        // Check single precision table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i)
            tot_str += ":"+col_sht.string(i);
        if (tot_str != sum_sht_str) {
            cout << endl << "TEST ERROR: GFitsTableShtCol::string()";
            cout << endl << "  Reference sum: " << sum_sht_str;
            cout << endl << "  Derived sum:   " << tot_str << endl;
            throw;
        }
        cout << ".";

        //
        // ===== S H O R T  1 0 =====
        //

        // Set short table
        GFitsTableShtCol col_sht10 = GFitsTableShtCol("SHORT10", nrows, nvec);
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j) {
                short val_sht  = short(100.0*cos(0.1*float(i))*
                                             cos(0.33*float(j)));
                col_sht10(i,j) = val_sht;
                sum_sht10     += val_sht;
                sum_sht10_int += int(val_sht);
                ostringstream s_value;
                s_value << val_sht;
                sum_sht10_str += ":"+s_value.str();
            }
        }
        cout << ".";

        // Check short table (operator access)
        tot_sht = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_sht += col_sht10(i,j);
        }
        if (tot_sht != sum_sht10) {
            cout << endl << "TEST ERROR: GFitsTableShtCol::operator() - 10";
            cout << endl << "  Reference sum: " << sum_sht10;
            cout << endl << "  Derived sum:   " << tot_sht << endl;
            throw;
        }
        cout << ".";

        // Check short table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_dbl += col_sht10.real(i,j);
        }
        if (!dequal(tot_dbl, sum_sht10)) {
            cout << endl << "TEST ERROR: GFitsTableShtCol::real() - 10";
            cout << endl << "  Reference sum: " << sum_sht10;
            cout << endl << "  Derived sum:   " << tot_dbl << endl;
            throw;
        }
        cout << ".";

        // Check short table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_int += col_sht10.integer(i,j);
        }
        if (tot_int != sum_sht10_int) {
            cout << endl << "TEST ERROR: GFitsTableShtCol::integer() - 10";
            cout << endl << "  Reference sum: " << sum_sht10_int;
            cout << endl << "  Derived sum:   " << tot_int << endl;
            throw;
        }
        cout << ".";

        // Check single precision table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_str += ":"+col_sht10.string(i,j);
        }
        if (tot_str != sum_sht10_str) {
            cout << endl << "TEST ERROR: GFitsTableShtCol::string() - 10";
            cout << endl << "  Reference sum: " << sum_sht10_str;
            cout << endl << "  Derived sum:   " << tot_str << endl;
            throw;
        }
        cout << ".";

        //
        // ===== W R I T E   T A B L E =====
        //

        // Build binary table
        GFitsBinTable table = GFitsBinTable(nrows);
        table.append_column(col_sht);
        table.insert_column(0, col_sht10);

        // Create HDU and append to FILE file
        GFitsHDU hdu(table);
        fits.append(hdu);

        // Save FITS file
        fits.save();
        cout << ".";

    }
    catch (exception &e) {
        cout << endl << "TEST ERROR: Unable to build tables." << endl;
        cout << e.what() << endl;
        throw;
    }
    cout << ".";


    //================================================================
    //================================================================
    //================================================================
    // Read tables
    //================================================================
    //================================================================
    //================================================================
    try {
        // Open FITS file
        GFits fits;
        fits.open(filename);

        // Set number of rows
        int nrows = 10;
        int nvec  = 10;

        // Initial control sums
        float       tot_flt = 0.0;
        short       tot_sht = 0;
        double      tot_dbl = 0.0;
        int         tot_int = 0;
        std::string tot_str;

        //
        // ===== S H O R T =====
        //

        // Get column
        GFitsTableShtCol* col_sht =
                    (GFitsTableShtCol*)fits.hdu(1)->column("SHORT");

        // Check short table (operator access)
        tot_sht = 0;
        for (int i = 0; i < nrows; ++i)
            tot_sht += (*col_sht)(i);
        if (tot_sht != sum_sht) {
            cout << endl << "TEST ERROR: GFitsTableShtCol::operator()";
            cout << endl << "  Reference sum: " << sum_sht;
            cout << endl << "  Derived sum:   " << tot_sht << endl;
            throw;
        }
        cout << ".";

        // Check short table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i)
            tot_dbl += col_sht->real(i);
        if (!dequal(tot_dbl, double(sum_sht))) {
            cout << endl << "TEST ERROR: GFitsTableShtCol::real()";
            cout << endl << "  Reference sum: " << sum_sht;
            cout << endl << "  Derived sum:   " << tot_dbl << endl;
            throw;
        }
        cout << ".";

        // Check short table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i)
            tot_int += col_sht->integer(i);
        if (tot_int != sum_sht_int) {
            cout << endl << "TEST ERROR: GFitsTableShtCol::integer()";
            cout << endl << "  Reference sum: " << sum_sht_int;
            cout << endl << "  Derived sum:   " << tot_int << endl;
            throw;
        }
        cout << ".";

        // Check single precision table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i)
            tot_str += ":"+col_sht->string(i);
        if (tot_str != sum_sht_str) {
            cout << endl << "TEST ERROR: GFitsTableShtCol::string()";
            cout << endl << "  Reference sum: " << sum_sht_str;
            cout << endl << "  Derived sum:   " << tot_str << endl;
            throw;
        }
        cout << ".";

        //
        // ===== S H O R T  1 0 =====
        //

        // Get column
        GFitsTableShtCol* col_sht10 =
                    (GFitsTableShtCol*)fits.hdu(1)->column("SHORT10");

        // Check short table (operator access)
        tot_sht = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_sht += (*col_sht10)(i,j);
        }
        if (tot_sht != sum_sht10) {
            cout << endl << "TEST ERROR: GFitsTableShtCol::operator() - 10";
            cout << endl << "  Reference sum: " << sum_sht10;
            cout << endl << "  Derived sum:   " << tot_sht << endl;
            throw;
        }
        cout << ".";

        // Check short table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_dbl += col_sht10->real(i,j);
        }
        if (!dequal(tot_dbl, sum_sht10)) {
            cout << endl << "TEST ERROR: GFitsTableShtCol::real() - 10";
            cout << endl << "  Reference sum: " << sum_sht10;
            cout << endl << "  Derived sum:   " << tot_dbl << endl;
            throw;
        }
        cout << ".";

        // Check short table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_int += col_sht10->integer(i,j);
        }
        if (tot_int != sum_sht10_int) {
            cout << endl << "TEST ERROR: GFitsTableShtCol::integer() - 10";
            cout << endl << "  Reference sum: " << sum_sht10_int;
            cout << endl << "  Derived sum:   " << tot_int << endl;
            throw;
        }
        cout << ".";

        // Check single precision table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_str += ":"+col_sht10->string(i,j);
        }
        if (tot_str != sum_sht10_str) {
            cout << endl << "TEST ERROR: GFitsTableShtCol::string() - 10";
            cout << endl << "  Reference sum: " << sum_sht10_str;
            cout << endl << "  Derived sum:   " << tot_str << endl;
            throw;
        }
        cout << ".";
    }
    catch (exception &e) {
        cout << endl << "TEST ERROR: Unable to build tables." << endl;
        cout << e.what() << endl;
        throw;
    }
    cout << ". ok." << endl;

}


/***************************************************************************
 *                       Test: FITS binary table test                      *
 ***************************************************************************/
void test_bintable_long(void)
{
    // Set filename
    std::string filename = "test_bintable_long.fits";

    // Dump header
    cout << "Test GFitsTableLngCol: ";
    
    // Remove FITS file
    std::string cmd = "rm -rf "+ filename;
    system(cmd.c_str());

    // Allocate reference sums
    long        sum_lng       = 0;
    long        sum_lng10     = 0;
    long        sum_lng_int   = 0;
    long        sum_lng10_int = 0;
    std::string sum_lng_str;
    std::string sum_lng10_str;

    // Build tables
    try {
        // Open FITS file
        GFits fits;
        fits.open(filename);

        // Set number of rows
        int nrows = 10;
        int nvec  = 10;

        // Initial control sums
        float       tot_flt = 0.0;
        double      tot_dbl = 0.0;
        int         tot_int = 0;
        long        tot_lng = 0;
        short       tot_sht = 0;
        std::string tot_str;

        //
        // ===== L O N G =====
        //

        // Set long table
        GFitsTableLngCol col_lng = GFitsTableLngCol("LONG", nrows);
        for (int i = 0; i < nrows; ++i) {
            long val_lng  = long(100000.0 * cos(0.1*float(i)));
            col_lng(i)    = val_lng;
            sum_lng      += val_lng;
            sum_lng_int  += int(val_lng);
            ostringstream s_value;
            s_value << val_lng;
            sum_lng_str += ":"+s_value.str();
        }
        cout << ".";

        // Check long table (operator access)
        tot_lng = 0;
        for (int i = 0; i < nrows; ++i)
            tot_lng += col_lng(i);
        if (tot_lng != sum_lng) {
            cout << endl << "TEST ERROR: GFitsTableLngCol::operator()";
            cout << endl << "  Reference sum: " << sum_lng;
            cout << endl << "  Derived sum:   " << tot_lng << endl;
            throw;
        }
        cout << ".";

        // Check long table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i)
            tot_dbl += col_lng.real(i);
        if (!dequal(tot_dbl, double(sum_lng))) {
            cout << endl << "TEST ERROR: GFitsTableLngCol::real()";
            cout << endl << "  Reference sum: " << sum_lng;
            cout << endl << "  Derived sum:   " << tot_dbl << endl;
            throw;
        }
        cout << ".";

        // Check long table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i)
            tot_int += col_lng.integer(i);
        if (tot_int != sum_lng_int) {
            cout << endl << "TEST ERROR: GFitsTableLngCol::integer()";
            cout << endl << "  Reference sum: " << sum_lng_int;
            cout << endl << "  Derived sum:   " << tot_int << endl;
            throw;
        }
        cout << ".";

        // Check long table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i)
            tot_str += ":"+col_lng.string(i);
        if (tot_str != sum_lng_str) {
            cout << endl << "TEST ERROR: GFitsTableLngCol::string()";
            cout << endl << "  Reference sum: " << sum_lng_str;
            cout << endl << "  Derived sum:   " << tot_str << endl;
            throw;
        }
        cout << ".";

        //
        // ===== L O N G  1 0 =====
        //

        // Set long table
        GFitsTableLngCol col_lng10 = GFitsTableLngCol("LONG10", nrows, nvec);
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j) {
                long val_lng   = long(1000.0*cos(0.1*float(i))*
                                             cos(0.33*float(j)));
                col_lng10(i,j) = val_lng;
                sum_lng10     += val_lng;
                sum_lng10_int += int(val_lng);
                ostringstream s_value;
                s_value << val_lng;
                sum_lng10_str += ":"+s_value.str();
            }
        }
        cout << ".";

        // Check long table (operator access)
        tot_lng = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_lng += col_lng10(i,j);
        }
        if (tot_lng != sum_lng10) {
            cout << endl << "TEST ERROR: GFitsTableLngCol::operator() - 10";
            cout << endl << "  Reference sum: " << sum_lng10;
            cout << endl << "  Derived sum:   " << tot_lng << endl;
            throw;
        }
        cout << ".";

        // Check long table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_dbl += col_lng10.real(i,j);
        }
        if (!dequal(tot_dbl, sum_lng10)) {
            cout << endl << "TEST ERROR: GFitsTableLngCol::real() - 10";
            cout << endl << "  Reference sum: " << sum_lng10;
            cout << endl << "  Derived sum:   " << tot_dbl << endl;
            throw;
        }
        cout << ".";

        // Check long table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_int += col_lng10.integer(i,j);
        }
        if (tot_int != sum_lng10_int) {
            cout << endl << "TEST ERROR: GFitsTableLngCol::integer() - 10";
            cout << endl << "  Reference sum: " << sum_lng10_int;
            cout << endl << "  Derived sum:   " << tot_int << endl;
            throw;
        }
        cout << ".";

        // Check single precision table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_str += ":"+col_lng10.string(i,j);
        }
        if (tot_str != sum_lng10_str) {
            cout << endl << "TEST ERROR: GFitsTableLngCol::string() - 10";
            cout << endl << "  Reference sum: " << sum_lng10_str;
            cout << endl << "  Derived sum:   " << tot_str << endl;
            throw;
        }
        cout << ".";

        //
        // ===== W R I T E   T A B L E =====
        //

        // Build binary table
        GFitsBinTable table = GFitsBinTable(nrows);
        table.append_column(col_lng);
        table.insert_column(99, col_lng10);

        // Create HDU and append to FILE file
        GFitsHDU hdu(table);
        fits.append(hdu);

        // Save FITS file
        fits.save();
        cout << ".";
    }
    catch (exception &e) {
        cout << endl << "TEST ERROR: Unable to build tables." << endl;
        cout << e.what() << endl;
        throw;
    }
    cout << ".";


    //================================================================
    //================================================================
    //================================================================
    // Read tables
    //================================================================
    //================================================================
    //================================================================
    try {
        // Open FITS file
        GFits fits;
        fits.open(filename);

        // Set number of rows
        int nrows = 10;
        int nvec  = 10;

        // Initial control sums
        float       tot_flt = 0.0;
        double      tot_dbl = 0.0;
        int         tot_int = 0;
        long        tot_lng = 0;
        std::string tot_str;

        //
        // ===== L O N G =====
        //

        // Get column
        GFitsTableLngCol* col_lng =
                    (GFitsTableLngCol*)fits.hdu(1)->column("LONG");

        // Check long table (operator access)
        tot_lng = 0;
        for (int i = 0; i < nrows; ++i)
            tot_lng += (*col_lng)(i);
        if (tot_lng != sum_lng) {
            cout << endl << "TEST ERROR: GFitsTableLngCol::operator()";
            cout << endl << "  Reference sum: " << sum_lng;
            cout << endl << "  Derived sum:   " << tot_lng << endl;
            throw;
        }
        cout << ".";

        // Check long table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i)
            tot_dbl += col_lng->real(i);
        if (!dequal(tot_dbl, double(sum_lng))) {
            cout << endl << "TEST ERROR: GFitsTableLngCol::real()";
            cout << endl << "  Reference sum: " << sum_lng;
            cout << endl << "  Derived sum:   " << tot_dbl << endl;
            throw;
        }
        cout << ".";

        // Check long table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i)
            tot_int += col_lng->integer(i);
        if (tot_int != sum_lng_int) {
            cout << endl << "TEST ERROR: GFitsTableLngCol::integer()";
            cout << endl << "  Reference sum: " << sum_lng_int;
            cout << endl << "  Derived sum:   " << tot_int << endl;
            throw;
        }
        cout << ".";

        // Check long table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i)
            tot_str += ":"+col_lng->string(i);
        if (tot_str != sum_lng_str) {
            cout << endl << "TEST ERROR: GFitsTableLngCol::string()";
            cout << endl << "  Reference sum: " << sum_lng_str;
            cout << endl << "  Derived sum:   " << tot_str << endl;
            throw;
        }
        cout << ".";

        //
        // ===== L O N G  1 0 =====
        //

        // Get column
        GFitsTableLngCol* col_lng10 =
                    (GFitsTableLngCol*)fits.hdu(1)->column("LONG10");

        // Check long table (operator access)
        tot_lng = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_lng += (*col_lng10)(i,j);
        }
        if (tot_lng != sum_lng10) {
            cout << endl << "TEST ERROR: GFitsTableLngCol::operator() - 10";
            cout << endl << "  Reference sum: " << sum_lng10;
            cout << endl << "  Derived sum:   " << tot_lng << endl;
            throw;
        }
        cout << ".";

        // Check long table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_dbl += col_lng10->real(i,j);
        }
        if (!dequal(tot_dbl, sum_lng10)) {
            cout << endl << "TEST ERROR: GFitsTableLngCol::real() - 10";
            cout << endl << "  Reference sum: " << sum_lng10;
            cout << endl << "  Derived sum:   " << tot_dbl << endl;
            throw;
        }
        cout << ".";

        // Check long table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_int += col_lng10->integer(i,j);
        }
        if (tot_int != sum_lng10_int) {
            cout << endl << "TEST ERROR: GFitsTableLngCol::integer() - 10";
            cout << endl << "  Reference sum: " << sum_lng10_int;
            cout << endl << "  Derived sum:   " << tot_int << endl;
            throw;
        }
        cout << ".";

        // Check single precision table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_str += ":"+col_lng10->string(i,j);
        }
        if (tot_str != sum_lng10_str) {
            cout << endl << "TEST ERROR: GFitsTableLngCol::string() - 10";
            cout << endl << "  Reference sum: " << sum_lng10_str;
            cout << endl << "  Derived sum:   " << tot_str << endl;
            throw;
        }
        cout << ".";
    }
    catch (exception &e) {
        cout << endl << "TEST ERROR: Unable to build tables." << endl;
        cout << e.what() << endl;
        throw;
    }
    cout << ". ok." << endl;

}


/***************************************************************************
 *                       Test: FITS binary table test                      *
 ***************************************************************************/
void test_bintable_string(void)
{
    // Set filename
    std::string filename = "test_bintable_string.fits";

    // Dump header
    cout << "Test GFitsTableStrCol: ";
    
    // Remove FITS file
    std::string cmd = "rm -rf "+ filename;
    system(cmd.c_str());

    // Allocate reference sums
    std::string sum_str;
    std::string sum_str10;
    double      sum_str_dbl   = 0;
    double      sum_str10_dbl = 0;
    double      sum_str_int   = 0;
    double      sum_str10_int = 0;
    
    // Build tables
    try {
        // Open FITS file
        GFits fits;
        fits.open(filename);

        // Set number of rows
        int nrows = 10;
        int nvec  = 10;

        // Initial control sums
        float       tot_flt = 0.0;
        double      tot_dbl = 0.0;
        int         tot_int = 0;
        long        tot_lng = 0;
        short       tot_sht = 0;
        std::string tot_str;

        //
        // ===== S T R I N G =====
        //

        // Set string table
        GFitsTableStrCol col_str = GFitsTableStrCol("STRING", nrows, 20);
        for (int i = 0; i < nrows; ++i) {
            double val_dbl = 100.0 * cos(0.1*double(i));
            ostringstream s_value;
            s_value << scientific << val_dbl;
            col_str(i)   = s_value.str();
            sum_str     += ":" + s_value.str();
            sum_str_dbl += val_dbl;
            sum_str_int += int(val_dbl);
        }
        cout << ".";

        // Check string table (operator access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i)
            tot_str += ":" + col_str(i);
        if (tot_str != sum_str) {
            cout << endl << "TEST ERROR: GFitsTableStrCol::operator()";
            cout << endl << "  Reference sum: " << sum_str;
            cout << endl << "  Derived sum:   " << tot_str << endl;
            throw;
        }
        cout << ".";

        // Check string table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i)
            tot_dbl += col_str.real(i);
        if (!equal(tot_dbl, sum_str_dbl, 1.0e-4)) {
            cout << endl << "TEST ERROR: GFitsTableStrCol::real()";
            cout << endl << "  Reference sum: " << sum_str_dbl;
            cout << endl << "  Derived sum:   " << tot_dbl << endl;
            throw;
        }
        cout << ".";

        // Check string table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i)
            tot_int += col_str.integer(i);
        if (tot_int != sum_str_int) {
            cout << endl << "TEST ERROR: GFitsTableStrCol::integer()";
            cout << endl << "  Reference sum: " << sum_str_int;
            cout << endl << "  Derived sum:   " << tot_int << endl;
            throw;
        }
        cout << ".";

        // Check string table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i)
            tot_str += ":"+col_str.string(i);
        if (tot_str != sum_str) {
            cout << endl << "TEST ERROR: GFitsTableStrCol::string()";
            cout << endl << "  Reference sum: " << sum_str;
            cout << endl << "  Derived sum:   " << tot_str << endl;
            throw;
        }
        cout << ".";

        //
        // ===== S T R I N G   1 0 =====
        //

        // Set string table
        GFitsTableStrCol col_str10 = GFitsTableStrCol("STRING10", nrows, 20, nvec);
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j) {
                double val_dbl = 100.0 * cos(0.1*double(i)) * cos(0.3*double(j));
                ostringstream s_value;
                s_value << scientific << val_dbl;
                col_str10(i,j)  = s_value.str();
                sum_str10      += ":" + s_value.str();
                sum_str10_dbl  += val_dbl;
                sum_str10_int  += int(val_dbl);
            }
        }
        cout << ".";

        // Check string table (operator access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_str += ":" + col_str10(i,j);
        }
        if (tot_str != sum_str10) {
            cout << endl << "TEST ERROR: GFitsTableStrCol::operator() - 10";
            cout << endl << "  Reference sum: " << sum_str10;
            cout << endl << "  Derived sum:   " << tot_str << endl;
            throw;
        }
        cout << ".";

        // Check string table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_dbl += col_str10.real(i,j);
        }
        if (!equal(tot_dbl, sum_str10_dbl, 1.0e-4)) {
            cout << endl << "TEST ERROR: GFitsTableStrCol::real() - 10";
            cout << endl << "  Reference sum: " << sum_str10_dbl;
            cout << endl << "  Derived sum:   " << tot_dbl << endl;
            throw;
        }
        cout << ".";

        // Check string table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_int += col_str10.integer(i,j);
        }
        if (tot_int != sum_str10_int) {
            cout << endl << "TEST ERROR: GFitsTableStrCol::integer() - 10";
            cout << endl << "  Reference sum: " << sum_str10_int;
            cout << endl << "  Derived sum:   " << tot_int << endl;
            throw;
        }
        cout << ".";

        // Check string table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_str += ":"+col_str10.string(i,j);
        }
        if (tot_str != sum_str10) {
            cout << endl << "TEST ERROR: GFitsTableStrCol::string() - 10";
            cout << endl << "  Reference sum: " << sum_str10;
            cout << endl << "  Derived sum:   " << tot_str << endl;
            throw;
        }
        cout << ".";

        //
        // ===== W R I T E   T A B L E =====
        //

        // Build binary table
        GFitsBinTable table = GFitsBinTable(nrows);
        table.append_column(col_str);
        table.append_column(col_str10);

        // Create HDU and append to FILE file
        GFitsHDU hdu(table);
        fits.append(hdu);

        // Save FITS file
        fits.save();
        cout << ".";
    }
    catch (exception &e) {
        cout << endl << "TEST ERROR: Unable to build tables." << endl;
        cout << e.what() << endl;
        throw;
    }
    cout << ".";


    //================================================================
    //================================================================
    //================================================================
    // Read tables
    //================================================================
    //================================================================
    //================================================================
    try {
        // Open FITS file
        GFits fits;
        fits.open(filename);

        // Set number of rows
        int nrows = 10;
        int nvec  = 10;

        // Initial control sums
        float       tot_flt = 0.0;
        double      tot_dbl = 0.0;
        int         tot_int = 0;
        std::string tot_str;

        //
        // ===== S T R I N G =====
        //

        // Get column
        GFitsTableStrCol* col_str =
                         (GFitsTableStrCol*)fits.hdu(1)->column("STRING");

        // Check string table (operator access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i)
            tot_str += ":" + (*col_str)(i);
        if (tot_str != sum_str) {
            cout << endl << "TEST ERROR: GFitsTableStrCol::operator()";
            cout << endl << "  Reference sum: " << sum_str;
            cout << endl << "  Derived sum:   " << tot_str << endl;
            throw;
        }
        cout << ".";

        // Check string table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i)
            tot_dbl += col_str->real(i);
        if (!equal(tot_dbl, sum_str_dbl, 1.0e-4)) {
            cout << endl << "TEST ERROR: GFitsTableStrCol::real()";
            cout << endl << "  Reference sum: " << sum_str_dbl;
            cout << endl << "  Derived sum:   " << tot_dbl << endl;
            throw;
        }
        cout << ".";

        // Check string table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i)
            tot_int += col_str->integer(i);
        if (tot_int != sum_str_int) {
            cout << endl << "TEST ERROR: GFitsTableStrCol::integer()";
            cout << endl << "  Reference sum: " << sum_str_int;
            cout << endl << "  Derived sum:   " << tot_int << endl;
            throw;
        }
        cout << ".";

        // Check string table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i)
            tot_str += ":"+col_str->string(i);
        if (tot_str != sum_str) {
            cout << endl << "TEST ERROR: GFitsTableStrCol::string()";
            cout << endl << "  Reference sum: " << sum_str;
            cout << endl << "  Derived sum:   " << tot_str << endl;
            throw;
        }
        cout << ".";

        //
        // ===== S T R I N G   1 0 =====
        //

        // Get column
        GFitsTableStrCol* col_str10 =
                         (GFitsTableStrCol*)fits.hdu(1)->column("STRING10");

        // Check string table (operator access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_str += ":" + (*col_str10)(i,j);
        }
        if (tot_str != sum_str10) {
            cout << endl << "TEST ERROR: GFitsTableStrCol::operator() - 10";
            cout << endl << "  Reference sum: " << sum_str10;
            cout << endl << "  Derived sum:   " << tot_str << endl;
            throw;
        }
        cout << ".";

        // Check string table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_dbl += col_str10->real(i,j);
        }
        if (!equal(tot_dbl, sum_str10_dbl, 1.0e-4)) {
            cout << endl << "TEST ERROR: GFitsTableStrCol::real() - 10";
            cout << endl << "  Reference sum: " << sum_str10_dbl;
            cout << endl << "  Derived sum:   " << tot_dbl << endl;
            throw;
        }
        cout << ".";

        // Check string table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_int += col_str10->integer(i,j);
        }
        if (tot_int != sum_str10_int) {
            cout << endl << "TEST ERROR: GFitsTableStrCol::integer() - 10";
            cout << endl << "  Reference sum: " << sum_str10_int;
            cout << endl << "  Derived sum:   " << tot_int << endl;
            throw;
        }
        cout << ".";

        // Check string table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_str += ":"+col_str10->string(i,j);
        }
        if (tot_str != sum_str10) {
            cout << endl << "TEST ERROR: GFitsTableStrCol::string() - 10";
            cout << endl << "  Reference sum: " << sum_str10;
            cout << endl << "  Derived sum:   " << tot_str << endl;
            throw;
        }
        cout << ".";

    }
    catch (exception &e) {
        cout << endl << "TEST ERROR: Unable to read tables." << endl;
        cout << e.what() << endl;
        throw;
    }
    cout << ". ok." << endl;

}


/***************************************************************************
 *                       Test: FITS binary table test                      *
 ***************************************************************************/
void test_bintable_logical(void)
{
    // Set filename
    std::string filename = "test_bintable_logical.fits";

    // Dump header
    cout << "Test GFitsTableLogCol: ";
    
    // Remove FITS file
    std::string cmd = "rm -rf "+ filename;
    system(cmd.c_str());

    // Allocate reference sums
    int         sum       = 0;
    int         sum10     = 0;
    int         sum_int   = 0;
    int         sum10_int = 0;
    std::string sum_str;
    std::string sum10_str;

    // Build tables
    try {
        // Open FITS file
        GFits fits;
        fits.open(filename);

        // Set number of rows
        int nrows = 10;
        int nvec  = 10;

        // Initial control sums
        float       tot_flt = 0.0;
        double      tot_dbl = 0.0;
        int         tot_int = 0;
        long        tot_lng = 0;
        std::string tot_str;

        //
        // ===== L O G I C A L =====
        //

        // Set table
        GFitsTableLogCol col = GFitsTableLogCol("LOGICAL", nrows);
        for (int i = 0; i < nrows; ++i) {
            char val  = (i % 2) ? 1 : 0;
            col(i)    = val;
            sum      += val;
            sum_int  += int(val);
            if (val)
                sum_str += ":T";
            else
                sum_str += ":F";
        }
        cout << ".";

        // Check table (operator access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i)
            tot_int += col(i);
        if (tot_int != sum) {
            cout << endl << "TEST ERROR: GFitsTableLogCol::operator()";
            cout << endl << "  Reference sum: " << sum;
            cout << endl << "  Derived sum:   " << tot_int << endl;
            throw;
        }
        cout << ".";

        // Check table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i)
            tot_dbl += col.real(i);
        if (!dequal(tot_dbl, double(sum))) {
            cout << endl << "TEST ERROR: GFitsTableLogCol::real()";
            cout << endl << "  Reference sum: " << sum;
            cout << endl << "  Derived sum:   " << tot_dbl << endl;
            throw;
        }
        cout << ".";

        // Check table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i)
            tot_int += col.integer(i);
        if (tot_int != sum_int) {
            cout << endl << "TEST ERROR: GFitsTableLogCol::integer()";
            cout << endl << "  Reference sum: " << sum_int;
            cout << endl << "  Derived sum:   " << tot_int << endl;
            throw;
        }
        cout << ".";

        // Check table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i)
            tot_str += ":"+col.string(i);
        if (tot_str != sum_str) {
            cout << endl << "TEST ERROR: GFitsTableLogCol::string()";
            cout << endl << "  Reference sum: " << sum_str;
            cout << endl << "  Derived sum:   " << tot_str << endl;
            throw;
        }
        cout << ".";

        //
        // ===== L O G I C A L  1 0 =====
        //

        // Set table
        GFitsTableLogCol col10 = GFitsTableLogCol("LOGICAL10", nrows, nvec);
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j) {
                char val  = (i % 2) * (j % 2);
                col10(i,j) = val;
                sum10     += val;
                sum10_int += int(val);
                if (val)
                    sum10_str += ":T";
                else
                    sum10_str += ":F";
            }
        }
        cout << ".";

        // Check table (operator access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_int += col10(i,j);
        }
        if (tot_int != sum10) {
            cout << endl << "TEST ERROR: GFitsTableLogCol::operator() - 10";
            cout << endl << "  Reference sum: " << sum10;
            cout << endl << "  Derived sum:   " << tot_int << endl;
            throw;
        }
        cout << ".";

        // Check table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_dbl += col10.real(i,j);
        }
        if (!dequal(tot_dbl, sum10)) {
            cout << endl << "TEST ERROR: GFitsTableLogCol::real() - 10";
            cout << endl << "  Reference sum: " << sum10;
            cout << endl << "  Derived sum:   " << tot_dbl << endl;
            throw;
        }
        cout << ".";

        // Check table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_int += col10.integer(i,j);
        }
        if (tot_int != sum10_int) {
            cout << endl << "TEST ERROR: GFitsTableLogCol::integer() - 10";
            cout << endl << "  Reference sum: " << sum10_int;
            cout << endl << "  Derived sum:   " << tot_int << endl;
            throw;
        }
        cout << ".";

        // Check table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_str += ":"+col10.string(i,j);
        }
        if (tot_str != sum10_str) {
            cout << endl << "TEST ERROR: GFitsTableLogCol::string() - 10";
            cout << endl << "  Reference sum: " << sum10_str;
            cout << endl << "  Derived sum:   " << tot_str << endl;
            throw;
        }
        cout << ".";

        //
        // ===== W R I T E   T A B L E =====
        //

        // Build binary table
        GFitsBinTable table = GFitsBinTable(nrows);
        table.append_column(col);
        table.insert_column(0, col10);

        // Create HDU and append to FILE file
        GFitsHDU hdu(table);
        fits.append(hdu);

        // Save FITS file
        fits.save();
        cout << ".";

    }
    catch (exception &e) {
        cout << endl << "TEST ERROR: Unable to build tables." << endl;
        cout << e.what() << endl;
        throw;
    }
    cout << ".";


    //================================================================
    //================================================================
    //================================================================
    // Read tables
    //================================================================
    //================================================================
    //================================================================
    try {
        // Open FITS file
        GFits fits;
        fits.open(filename);

        // Set number of rows
        int nrows = 10;
        int nvec  = 10;

        // Initial control sums
        float       tot_flt = 0.0;
        short       tot_sht = 0;
        double      tot_dbl = 0.0;
        int         tot_int = 0;
        std::string tot_str;

        //
        // ===== L O G I C A L =====
        //

        // Get column
        GFitsTableLogCol* col =
                    (GFitsTableLogCol*)fits.hdu(1)->column("LOGICAL");

        // Check table (operator access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i)
            tot_int += (*col)(i);
        if (tot_int != sum) {
            cout << endl << "TEST ERROR: GFitsTableLogCol::operator()";
            cout << endl << "  Reference sum: " << sum;
            cout << endl << "  Derived sum:   " << tot_int << endl;
            throw;
        }
        cout << ".";

        // Check table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i)
            tot_dbl += col->real(i);
        if (!dequal(tot_dbl, double(sum))) {
            cout << endl << "TEST ERROR: GFitsTableLogCol::real()";
            cout << endl << "  Reference sum: " << sum;
            cout << endl << "  Derived sum:   " << tot_dbl << endl;
            throw;
        }
        cout << ".";

        // Check table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i)
            tot_int += col->integer(i);
        if (tot_int != sum_int) {
            cout << endl << "TEST ERROR: GFitsTableLogCol::integer()";
            cout << endl << "  Reference sum: " << sum_int;
            cout << endl << "  Derived sum:   " << tot_int << endl;
            throw;
        }
        cout << ".";

        // Check single precision table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i)
            tot_str += ":"+col->string(i);
        if (tot_str != sum_str) {
            cout << endl << "TEST ERROR: GFitsTableLogCol::string()";
            cout << endl << "  Reference sum: " << sum_str;
            cout << endl << "  Derived sum:   " << tot_str << endl;
            throw;
        }
        cout << ".";

        //
        // ===== L O G I C A L  1 0 =====
        //

        // Get column
        GFitsTableLogCol* col10 =
                    (GFitsTableLogCol*)fits.hdu(1)->column("LOGICAL10");

        // Check table (operator access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_int += (*col10)(i,j);
        }
        if (tot_int != sum10) {
            cout << endl << "TEST ERROR: GFitsTableLogCol::operator() - 10";
            cout << endl << "  Reference sum: " << sum10;
            cout << endl << "  Derived sum:   " << tot_int << endl;
            throw;
        }
        cout << ".";

        // Check table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_dbl += col10->real(i,j);
        }
        if (!dequal(tot_dbl, sum10)) {
            cout << endl << "TEST ERROR: GFitsTableLogCol::real() - 10";
            cout << endl << "  Reference sum: " << sum10;
            cout << endl << "  Derived sum:   " << tot_dbl << endl;
            throw;
        }
        cout << ".";

        // Check table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_int += col10->integer(i,j);
        }
        if (tot_int != sum10_int) {
            cout << endl << "TEST ERROR: GFitsTableLogCol::integer() - 10";
            cout << endl << "  Reference sum: " << sum10_int;
            cout << endl << "  Derived sum:   " << tot_int << endl;
            throw;
        }
        cout << ".";

        // Check single precision table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_str += ":"+col10->string(i,j);
        }
        if (tot_str != sum10_str) {
            cout << endl << "TEST ERROR: GFitsTableLogCol::string() - 10";
            cout << endl << "  Reference sum: " << sum10_str;
            cout << endl << "  Derived sum:   " << tot_str << endl;
            throw;
        }
        cout << ".";
    }
    catch (exception &e) {
        cout << endl << "TEST ERROR: Unable to build tables." << endl;
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
    test_create();
    test_image_double();
    test_bintable_double();
    test_bintable_float();
    test_bintable_short();
    test_bintable_long();
    test_bintable_string();
    test_bintable_logical();
    
    // Return
    return 0;
}
