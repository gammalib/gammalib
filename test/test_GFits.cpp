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
#include <iostream>                           // std::cout, cerr
#include <stdexcept>                          // std::exception
#include <stdlib.h>
#include "GammaLib.hpp"
#include "GTools.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Globals ____________________________________________________________ */


/***************************************************************************
 * @brief Verification routines
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
 * @brief Test FITS file creation
 ***************************************************************************/
void test_create(void)
{
    // Dump header
    std::cout << "Test GFits: Create FITS file: ";

    // Remove FITS file
    system("rm -rf test_empty.fits");
    system("rm -rf test_empty_image.fits");
    system("rm -rf test.fits");
    system("rm -rf test_create_bintable.fits");

    // Create empty FITS file
    try {
        GFits fits("test_empty.fits");
        fits.save();
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to create empty FITS file."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Create FITS file with empty double precision image
    try {
        GFits fits;
        fits.open("test_empty_image.fits");
        GFitsImageDouble image;
        fits.append(&image);
        fits.save();
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to create FITS file with empty image."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Attach double precision image
    double sum = 0.0;
    try {
        GFits fits;
        fits.open("test_empty_image.fits");
        int naxis       = 2;
        int nx          = 10;
        int ny          = 20;
        int naxes[]     = {nx,ny};
        GFitsImageDouble image(naxis, naxes);
        for (int ix = 0; ix < nx; ++ix) {
            for (int iy = 0; iy < ny; ++iy) {
                image(ix,iy) = 0.01*ix + 0.01*iy;
                sum += image(ix,iy);
            }
        }
        fits.append(&image);
        fits.saveto("test.fits");
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to create FITS file."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Re-open double precision image
    try {
        GFits fits;
        fits.open("test.fits");
        GFitsHDU*         hdu   = fits.hdu(0);
        GFitsImageDouble* image = (GFitsImageDouble*)fits.hdu(1);
        int nx = image->naxes(0);
        int ny = image->naxes(1);
        double total = 0.0;
        for (int ix = 0; ix < nx; ++ix) {
            for (int iy = 0; iy < ny; ++iy) {
                total += (*image)(ix,iy);
            }
        }
        if (!dequal(total, sum)) {
            std::cout << std::endl
                      << "TEST ERROR: Bad values in loaded image."
                      << std::endl;
            throw;
        }
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to re-open FITS file."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Attach binary table (save variant)
    try {
        // Re-open FITS file
        GFits fits;
        fits.open("test.fits");

        // Create binary Table with 10 rows
        int nrows = 10;
        GFitsBinTable       table  = GFitsBinTable(nrows);
        GFitsTableDoubleCol first  = GFitsTableDoubleCol("First", nrows);
        GFitsTableDoubleCol second = GFitsTableDoubleCol("Second", nrows);
        first(0)  =  1.0;
        first(1)  =  2.0;
        first(2)  =  3.0;
        second(0) = -99.0;
        //std::cout << std::endl << first << std::endl;
        //std::cout << second << std::endl;
        table.append_column(first);
        table.append_column(second);

        // Append table to FILE file
        fits.append(&table);

        // Save FITS file
        fits.save();
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to attach binary table."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Crate binary table (saveto variant)
    try {
        // Allocate FITS file
        GFits fits;

        // Create binary Table with 10 rows
        int nrows = 10;
        GFitsBinTable       table  = GFitsBinTable(nrows);
        GFitsTableDoubleCol first  = GFitsTableDoubleCol("First", nrows);
        GFitsTableDoubleCol second = GFitsTableDoubleCol("Second", nrows);
        first(0)  =  1.0;
        first(1)  =  2.0;
        first(2)  =  3.0;
        second(0) = -99.0;
        table.append_column(first);
        table.append_column(second);

        // Append table to FILE file
        fits.append(&table);

        // Save FITS file
        fits.saveto("test_create_bintable.fits");
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to attach binary table."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

    // Final ok
    std::cout << ". ok." << std::endl;

    // Return
    return;
}


/***************************************************************************
 * @brief Test GFitsImageDouble class
 ***************************************************************************/
void test_image_double(void)
{
    // Dump header
    std::cout << "Test GFitsImageDouble: ";

    // Test 1D Image
    try {
        int naxis   = 1;
        int nx      = 10;
        int naxes[] = {nx};
        GFitsImageDouble image(naxis, naxes);
        std::cout << ".";

        // Fill and read image
        try {
            for (int ix = 0; ix < nx; ++ix)
                image(ix) = double(ix);
            double sum = 0.0;
            for (int ix = 0; ix < nx; ++ix)
                sum += image(ix);
            if (!dequal(sum, 0.5*double(nx-1)*double(nx))) {
                std::cout << std::endl
                          << "TEST ERROR: Bad fill/read of 1D image."
                          << std::endl;
                throw;
            }
        }
        catch (std::exception &e) {
            std::cout << std::endl
                      << "TEST ERROR: Unable to fill/read 1D image."
                      << std::endl;
            std::cout << e.what() << std::endl;
            throw;
        }
        std::cout << ".";

        // Bad 2D operator
        try {
            double result = image.at(0,1);
        }
        catch (GException::fits_wrong_image_operator &e) {
            std::cout << ".";
        }
        catch (std::exception &e) {
            std::cout << std::endl
                      << "TEST ERROR: Wrong operactor access error."
                      << std::endl;
            std::cout << e.what() << std::endl;
            throw;
        }

        // Bad 3D operator
        try {
            double result = image.at(0,1,2);
        }
        catch (GException::fits_wrong_image_operator &e) {
            std::cout << ".";
        }
        catch (std::exception &e) {
            std::cout << std::endl
                      << "TEST ERROR: Wrong operactor access error."
                      << std::endl;
            std::cout << e.what() << std::endl;
            throw;
        }

        // Bad 4D operator
        try {
            double result = image.at(0,1,2,3);
        }
        catch (GException::fits_wrong_image_operator &e) {
            std::cout << ".";
        }
        catch (std::exception &e) {
            std::cout << std::endl
                      << "TEST ERROR: Wrong operactor access error."
                      << std::endl;
            std::cout << e.what() << std::endl;
            throw;
        }

    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to test GFitsImageDouble classes."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ". ok." << std::endl;
}


/***************************************************************************
 * @brief Test double precision FITS binary table
 ***************************************************************************/
void test_bintable_double(void)
{
    // Set filename
    std::string filename = "test_bintable_double.fits";

    // Dump header
    std::cout << "Test GFitsTableDoubleCol: ";
    
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
        GFitsTableDoubleCol col_dbl = GFitsTableDoubleCol("DOUBLE", nrows);
        for (int i = 0; i < nrows; ++i) {
            double val_dbl = cos(double(i));
            col_dbl(i)   = val_dbl;
            sum_dbl     += val_dbl;
            sum_dbl_int += int(val_dbl);
            sum_dbl_str += ":"+str(val_dbl);
        }
        std::cout << ".";

        // Check double precision table (operator access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i)
            tot_dbl += col_dbl(i);
        if (!dequal(tot_dbl, sum_dbl)) {
            std::cout << std::endl << "TEST ERROR: GFitsTableDoubleCol::operator()";
            std::cout << std::endl << "  Reference sum: " << sum_dbl;
            std::cout << std::endl << "  Derived sum:   " << tot_dbl << std::endl;
            throw;
        }
        std::cout << ".";

        // Check double precision table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i)
            tot_dbl += col_dbl.real(i);
        if (!dequal(tot_dbl, sum_dbl)) {
            std::cout << std::endl << "TEST ERROR: GFitsTableDoubleCol::real()";
            std::cout << std::endl << "  Reference sum: " << sum_dbl;
            std::cout << std::endl << "  Derived sum:   " << tot_dbl << std::endl;
            throw;
        }
        std::cout << ".";

        // Check double precision table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i)
            tot_int += col_dbl.integer(i);
        if (tot_int != sum_dbl_int) {
            std::cout << std::endl << "TEST ERROR: GFitsTableDoubleCol::integer()";
            std::cout << std::endl << "  Reference sum: " << sum_dbl_int;
            std::cout << std::endl << "  Derived sum:   " << tot_int << std::endl;
            throw;
        }
        std::cout << ".";

        // Check double precision table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i)
            tot_str += ":"+col_dbl.string(i);
        if (tot_str != sum_dbl_str) {
            std::cout << std::endl << "TEST ERROR: GFitsTableDoubleCol::string()";
            std::cout << std::endl << "  Reference sum: " << sum_dbl_str;
            std::cout << std::endl << "  Derived sum:   " << tot_str << std::endl;
            throw;
        }
        std::cout << ".";

        //
        // ===== D O U B L E 1 0 =====
        //

        // Set double precision table
        GFitsTableDoubleCol col_dbl10 = GFitsTableDoubleCol("DOUBLE10", nrows, nvec);
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j) {
                double val_dbl = cos(double(i))*cos(0.33*double(j));
                col_dbl10(i,j) = val_dbl;
                sum_dbl10     += val_dbl;
                sum_dbl10_int += int(val_dbl);
                sum_dbl10_str += ":"+str(val_dbl);
            }
        }
        std::cout << ".";

        // Check double precision table (operator access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_dbl += col_dbl10(i,j);
        }
        if (!dequal(tot_dbl, sum_dbl10)) {
            std::cout << std::endl << "TEST ERROR: GFitsTableDoubleCol::operator() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_dbl10;
            std::cout << std::endl << "  Derived sum:   " << tot_dbl << std::endl;
            throw;
        }
        std::cout << ".";

        // Check double precision table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_dbl += col_dbl10.real(i,j);
        }
        if (!dequal(tot_dbl, sum_dbl10)) {
            std::cout << std::endl << "TEST ERROR: GFitsTableDoubleCol::real() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_dbl10;
            std::cout << std::endl << "  Derived sum:   " << tot_dbl << std::endl;
            throw;
        }
        std::cout << ".";

        // Check double precision table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_int += col_dbl10.integer(i,j);
        }
        if (tot_int != sum_dbl10_int) {
            std::cout << std::endl << "TEST ERROR: GFitsTableDoubleCol::integer() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_dbl10_int;
            std::cout << std::endl << "  Derived sum:   " << tot_int << std::endl;
            throw;
        }
        std::cout << ".";

        // Check double precision table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_str += ":"+col_dbl10.string(i,j);
        }
        if (tot_str != sum_dbl10_str) {
            std::cout << std::endl << "TEST ERROR: GFitsTableDoubleCol::string() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_dbl10_str;
            std::cout << std::endl << "  Derived sum:   " << tot_str << std::endl;
            throw;
        }
        std::cout << ".";

        //
        // ===== W R I T E   T A B L E =====
        //

        // Build binary table
        GFitsBinTable table = GFitsBinTable(nrows);
        table.append_column(col_dbl);
        table.insert_column(1, col_dbl10);

        // Append to FILE file
        fits.append(&table);

        // Save FITS file
        fits.save();
        std::cout << ".";
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to build GFitsTableDoubleCol tables."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";


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
        GFitsTableDoubleCol* col_dbl =
            (GFitsTableDoubleCol*)((GFitsTable*)fits.hdu(1))->column("DOUBLE");

        // Check double precision table (operator access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i)
            tot_dbl += (*col_dbl)(i);
        if (!dequal(tot_dbl, sum_dbl)) {
            std::cout << std::endl << "TEST ERROR: GFitsTableDoubleCol::operator()";
            std::cout << std::endl << "  Reference sum: " << sum_dbl;
            std::cout << std::endl << "  Derived sum:   " << tot_dbl << std::endl;
            throw;
        }
        std::cout << ".";

        // Check double precision table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i)
            tot_dbl += (*col_dbl).real(i);
        if (!dequal(tot_dbl, sum_dbl)) {
            std::cout << std::endl << "TEST ERROR: GFitsTableDoubleCol::real()";
            std::cout << std::endl << "  Reference sum: " << sum_dbl;
            std::cout << std::endl << "  Derived sum:   " << tot_dbl << std::endl;
            throw;
        }
        std::cout << ".";

        // Check double precision table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i)
            tot_int += (*col_dbl).integer(i);
        if (tot_int != sum_dbl_int) {
            std::cout << std::endl << "TEST ERROR: GFitsTableDoubleCol::integer()";
            std::cout << std::endl << "  Reference sum: " << sum_dbl_int;
            std::cout << std::endl << "  Derived sum:   " << tot_int << std::endl;
            throw;
        }
        std::cout << ".";

        // Check double precision table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i)
            tot_str += ":"+ (*col_dbl).string(i);
        if (tot_str != sum_dbl_str) {
            std::cout << std::endl << "TEST ERROR: GFitsTableDoubleCol::string()";
            std::cout << std::endl << "  Reference sum: " << sum_dbl_str;
            std::cout << std::endl << "  Derived sum:   " << tot_str << std::endl;
            throw;
        }
        std::cout << ".";

        //
        // ===== D O U B L E 1 0 =====
        //

        // Get column
        GFitsTableDoubleCol* col_dbl10 =
            (GFitsTableDoubleCol*)((GFitsTable*)fits.hdu(1))->column("DOUBLE10");

        // Check double precision table (operator access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_dbl += (*col_dbl10)(i,j);
        }
        if (!dequal(tot_dbl, sum_dbl10)) {
            std::cout << std::endl << "TEST ERROR: GFitsTableDoubleCol::operator() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_dbl10;
            std::cout << std::endl << "  Derived sum:   " << tot_dbl << std::endl;
            throw;
        }
        std::cout << ".";

        // Check double precision table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_dbl += col_dbl10->real(i,j);
        }
        if (!dequal(tot_dbl, sum_dbl10)) {
            std::cout << std::endl << "TEST ERROR: GFitsTableDoubleCol::real() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_dbl10;
            std::cout << std::endl << "  Derived sum:   " << tot_dbl << std::endl;
            throw;
        }
        std::cout << ".";

        // Check double precision table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_int += col_dbl10->integer(i,j);
        }
        if (tot_int != sum_dbl10_int) {
            std::cout << std::endl << "TEST ERROR: GFitsTableDoubleCol::integer() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_dbl10_int;
            std::cout << std::endl << "  Derived sum:   " << tot_int << std::endl;
            throw;
        }
        std::cout << ".";

        // Check double precision table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_str += ":"+col_dbl10->string(i,j);
        }
        if (tot_str != sum_dbl10_str) {
            std::cout << std::endl << "TEST ERROR: GFitsTableDoubleCol::string() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_dbl10_str;
            std::cout << std::endl << "  Derived sum:   " << tot_str << std::endl;
            throw;
        }
        std::cout << ".";
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to read GFitsTableDoubleCol tables."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ". ok." << std::endl;

}


/***************************************************************************
 * @brief Test single precision FITS binary table
 ***************************************************************************/
void test_bintable_float(void)
{
    // Set filename
    std::string filename = "test_bintable_float.fits";

    // Dump header
    std::cout << "Test GFitsTableFloatCol: ";
    
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
        GFitsTableFloatCol col_flt = GFitsTableFloatCol("FLOAT", nrows);
        for (int i = 0; i < nrows; ++i) {
            float val_flt = cos(0.1*float(i));
            col_flt(i)   = val_flt;
            sum_flt     += val_flt;
            sum_flt_int += int(val_flt);
            sum_flt_str += ":"+str(val_flt);
        }
        std::cout << ".";

        // Check single precision table (operator access)
        tot_flt = 0.0;
        for (int i = 0; i < nrows; ++i)
            tot_flt += col_flt(i);
        if (!fequal(tot_flt, sum_flt)) {
            std::cout << std::endl << "TEST ERROR: GFitsTableFloatCol::operator()";
            std::cout << std::endl << "  Reference sum: " << sum_flt;
            std::cout << std::endl << "  Derived sum:   " << tot_flt << std::endl;
            throw;
        }
        std::cout << ".";

        // Check single precision table (real access)
        tot_flt = 0.0;
        for (int i = 0; i < nrows; ++i)
            tot_flt += col_flt.real(i);
        if (!fequal(tot_flt, sum_flt)) {
            std::cout << std::endl << "TEST ERROR: GFitsTableFloatCol::real()";
            std::cout << std::endl << "  Reference sum: " << sum_flt;
            std::cout << std::endl << "  Derived sum:   " << tot_flt << std::endl;
            throw;
        }
        std::cout << ".";

        // Check single precision table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i)
            tot_int += col_flt.integer(i);
        if (tot_int != sum_flt_int) {
            std::cout << std::endl << "TEST ERROR: GFitsTableFloatCol::integer()";
            std::cout << std::endl << "  Reference sum: " << sum_flt_int;
            std::cout << std::endl << "  Derived sum:   " << tot_int << std::endl;
            throw;
        }
        std::cout << ".";

        // Check single precision table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i)
            tot_str += ":"+col_flt.string(i);
        if (tot_str != sum_flt_str) {
            std::cout << std::endl << "TEST ERROR: GFitsTableFloatCol::string()";
            std::cout << std::endl << "  Reference sum: " << sum_flt_str;
            std::cout << std::endl << "  Derived sum:   " << tot_str << std::endl;
            throw;
        }
        std::cout << ".";

        //
        // ===== F L O A T  1 0 =====
        //

        // Set single precision table
        GFitsTableFloatCol col_flt10 = GFitsTableFloatCol("FLOAT10", nrows, nvec);
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j) {
                float val_flt  = cos(0.1*float(i))*cos(0.33*float(j));
                col_flt10(i,j) = val_flt;
                sum_flt10     += val_flt;
                sum_flt10_int += int(val_flt);
                sum_flt10_str += ":"+str(val_flt);
            }
        }
        std::cout << ".";

        // Check single precision table (operator access)
        tot_flt = 0.0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_flt += col_flt10(i,j);
        }
        if (!fequal(tot_flt, sum_flt10)) {
            std::cout << std::endl << "TEST ERROR: GFitsTableFloatCol::operator() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_flt10;
            std::cout << std::endl << "  Derived sum:   " << tot_flt << std::endl;
            throw;
        }
        std::cout << ".";

        // Check single precision table (real access)
        tot_flt = 0.0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_flt += col_flt10.real(i,j);
        }
        if (!fequal(tot_flt, sum_flt10)) {
            std::cout << std::endl << "TEST ERROR: GFitsTableFloatCol::real() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_flt10;
            std::cout << std::endl << "  Derived sum:   " << tot_flt << std::endl;
            throw;
        }
        std::cout << ".";

        // Check double precision table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_int += col_flt10.integer(i,j);
        }
        if (tot_int != sum_flt10_int) {
            std::cout << std::endl << "TEST ERROR: GFitsTableFloatCol::integer() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_flt10_int;
            std::cout << std::endl << "  Derived sum:   " << tot_int << std::endl;
            throw;
        }
        std::cout << ".";

        // Check single precision table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_str += ":"+col_flt10.string(i,j);
        }
        if (tot_str != sum_flt10_str) {
            std::cout << std::endl << "TEST ERROR: GFitsTableFloatCol::string() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_flt10_str;
            std::cout << std::endl << "  Derived sum:   " << tot_str << std::endl;
            throw;
        }
        std::cout << ".";

        //
        // ===== W R I T E   T A B L E =====
        //

        // Build binary table
        GFitsBinTable table = GFitsBinTable(nrows);
        table.append_column(col_flt10);
        table.insert_column(0, col_flt);

        // Create HDU and append to FILE file
        fits.append(&table);

        // Save FITS file
        fits.save();
        std::cout << ".";
    }
    catch (std::exception &e) {
        std::cout << std::endl << "TEST ERROR: Unable to test floating point column in"
                        " binary tables." << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";


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
        GFitsTableFloatCol* col_flt =
            (GFitsTableFloatCol*)((GFitsTable*)fits.hdu(1))->column("FLOAT");

        // Check single precision table (operator access)
        tot_flt = 0.0;
        for (int i = 0; i < nrows; ++i)
            tot_flt += (*col_flt)(i);
        if (!fequal(tot_flt, sum_flt)) {
            std::cout << std::endl << "TEST ERROR: GFitsTableFloatCol::operator()";
            std::cout << std::endl << "  Reference sum: " << sum_flt;
            std::cout << std::endl << "  Derived sum:   " << tot_flt << std::endl;
            throw;
        }
        std::cout << ".";

        // Check single precision table (real access)
        tot_flt = 0.0;
        for (int i = 0; i < nrows; ++i)
            tot_flt += col_flt->real(i);
        if (!fequal(tot_flt, sum_flt)) {
            std::cout << std::endl << "TEST ERROR: GFitsTableFloatCol::real()";
            std::cout << std::endl << "  Reference sum: " << sum_flt;
            std::cout << std::endl << "  Derived sum:   " << tot_flt << std::endl;
            throw;
        }
        std::cout << ".";

        // Check single precision table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i)
            tot_int += col_flt->integer(i);
        if (tot_int != sum_flt_int) {
            std::cout << std::endl << "TEST ERROR: GFitsTableFloatCol::integer()";
            std::cout << std::endl << "  Reference sum: " << sum_flt_int;
            std::cout << std::endl << "  Derived sum:   " << tot_int << std::endl;
            throw;
        }
        std::cout << ".";

        // Check single precision table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i)
            tot_str += ":"+col_flt->string(i);
        if (tot_str != sum_flt_str) {
            std::cout << std::endl << "TEST ERROR: GFitsTableFloatCol::string()";
            std::cout << std::endl << "  Reference sum: " << sum_flt_str;
            std::cout << std::endl << "  Derived sum:   " << tot_str << std::endl;
            throw;
        }
        std::cout << ".";

        //
        // ===== F L O A T  1 0 =====
        //

        // Get column
        GFitsTableFloatCol* col_flt10 =
            (GFitsTableFloatCol*)((GFitsTable*)fits.hdu(1))->column("FLOAT10");

        // Check single precision table (operator access)
        tot_flt = 0.0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_flt += (*col_flt10)(i,j);
        }
        if (!fequal(tot_flt, sum_flt10)) {
            std::cout << std::endl << "TEST ERROR: GFitsTableFloatCol::operator() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_flt10;
            std::cout << std::endl << "  Derived sum:   " << tot_flt << std::endl;
            throw;
        }
        std::cout << ".";

        // Check single precision table (real access)
        tot_flt = 0.0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_flt += col_flt10->real(i,j);
        }
        if (!fequal(tot_flt, sum_flt10)) {
            std::cout << std::endl << "TEST ERROR: GFitsTableFloatCol::real() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_flt10;
            std::cout << std::endl << "  Derived sum:   " << tot_flt << std::endl;
            throw;
        }
        std::cout << ".";

        // Check double precision table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_int += col_flt10->integer(i,j);
        }
        if (tot_int != sum_flt10_int) {
            std::cout << std::endl << "TEST ERROR: GFitsTableFloatCol::integer() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_flt10_int;
            std::cout << std::endl << "  Derived sum:   " << tot_int << std::endl;
            throw;
        }
        std::cout << ".";

        // Check single precision table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_str += ":"+col_flt10->string(i,j);
        }
        if (tot_str != sum_flt10_str) {
            std::cout << std::endl << "TEST ERROR: GFitsTableFloatCol::string() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_flt10_str;
            std::cout << std::endl << "  Derived sum:   " << tot_str << std::endl;
            throw;
        }
        std::cout << ".";

    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to read GFitsTableFloatCol tables."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ". ok." << std::endl;

}


/***************************************************************************
 * @brief Test short FITS binary table
 ***************************************************************************/
void test_bintable_short(void)
{
    // Set filename
    std::string filename = "test_bintable_short.fits";

    // Dump header
    std::cout << "Test GFitsTableShortCol: ";
    
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
        GFitsTableShortCol col_sht = GFitsTableShortCol("SHORT", nrows);
        for (int i = 0; i < nrows; ++i) {
            short val_sht = short(1000.0 * cos(0.1*float(i)));
            col_sht(i)    = val_sht;
            sum_sht      += val_sht;
            sum_sht_int  += int(val_sht);
            sum_sht_str += ":"+str(val_sht);
        }
        std::cout << ".";

        // Check short table (operator access)
        tot_sht = 0;
        for (int i = 0; i < nrows; ++i)
            tot_sht += col_sht(i);
        if (tot_sht != sum_sht) {
            std::cout << std::endl << "TEST ERROR: GFitsTableShortCol::operator()";
            std::cout << std::endl << "  Reference sum: " << sum_sht;
            std::cout << std::endl << "  Derived sum:   " << tot_sht << std::endl;
            throw;
        }
        std::cout << ".";

        // Check short table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i)
            tot_dbl += col_sht.real(i);
        if (!dequal(tot_dbl, double(sum_sht))) {
            std::cout << std::endl << "TEST ERROR: GFitsTableShortCol::real()";
            std::cout << std::endl << "  Reference sum: " << sum_sht;
            std::cout << std::endl << "  Derived sum:   " << tot_dbl << std::endl;
            throw;
        }
        std::cout << ".";

        // Check short table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i)
            tot_int += col_sht.integer(i);
        if (tot_int != sum_sht_int) {
            std::cout << std::endl << "TEST ERROR: GFitsTableShortCol::integer()";
            std::cout << std::endl << "  Reference sum: " << sum_sht_int;
            std::cout << std::endl << "  Derived sum:   " << tot_int << std::endl;
            throw;
        }
        std::cout << ".";

        // Check single precision table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i)
            tot_str += ":"+col_sht.string(i);
        if (tot_str != sum_sht_str) {
            std::cout << std::endl << "TEST ERROR: GFitsTableShortCol::string()";
            std::cout << std::endl << "  Reference sum: " << sum_sht_str;
            std::cout << std::endl << "  Derived sum:   " << tot_str << std::endl;
            throw;
        }
        std::cout << ".";

        //
        // ===== S H O R T  1 0 =====
        //

        // Set short table
        GFitsTableShortCol col_sht10 = GFitsTableShortCol("SHORT10", nrows, nvec);
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j) {
                short val_sht  = short(100.0*cos(0.1*float(i))*
                                             cos(0.33*float(j)));
                col_sht10(i,j) = val_sht;
                sum_sht10     += val_sht;
                sum_sht10_int += int(val_sht);
                sum_sht10_str += ":"+str(val_sht);
            }
        }
        std::cout << ".";

        // Check short table (operator access)
        tot_sht = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_sht += col_sht10(i,j);
        }
        if (tot_sht != sum_sht10) {
            std::cout << std::endl << "TEST ERROR: GFitsTableShortCol::operator() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_sht10;
            std::cout << std::endl << "  Derived sum:   " << tot_sht << std::endl;
            throw;
        }
        std::cout << ".";

        // Check short table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_dbl += col_sht10.real(i,j);
        }
        if (!dequal(tot_dbl, sum_sht10)) {
            std::cout << std::endl << "TEST ERROR: GFitsTableShortCol::real() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_sht10;
            std::cout << std::endl << "  Derived sum:   " << tot_dbl << std::endl;
            throw;
        }
        std::cout << ".";

        // Check short table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_int += col_sht10.integer(i,j);
        }
        if (tot_int != sum_sht10_int) {
            std::cout << std::endl << "TEST ERROR: GFitsTableShortCol::integer() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_sht10_int;
            std::cout << std::endl << "  Derived sum:   " << tot_int << std::endl;
            throw;
        }
        std::cout << ".";

        // Check single precision table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_str += ":"+col_sht10.string(i,j);
        }
        if (tot_str != sum_sht10_str) {
            std::cout << std::endl << "TEST ERROR: GFitsTableShortCol::string() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_sht10_str;
            std::cout << std::endl << "  Derived sum:   " << tot_str << std::endl;
            throw;
        }
        std::cout << ".";

        //
        // ===== W R I T E   T A B L E =====
        //

        // Build binary table
        GFitsBinTable table = GFitsBinTable(nrows);
        table.append_column(col_sht);
        table.insert_column(0, col_sht10);

        // Append to file
        fits.append(&table);

        // Save FITS file
        fits.save();
        std::cout << ".";

    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to build GFitsTableShortCol tables."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";


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
        GFitsTableShortCol* col_sht =
            (GFitsTableShortCol*)((GFitsTable*)fits.hdu(1))->column("SHORT");

        // Check short table (operator access)
        tot_sht = 0;
        for (int i = 0; i < nrows; ++i)
            tot_sht += (*col_sht)(i);
        if (tot_sht != sum_sht) {
            std::cout << std::endl << "TEST ERROR: GFitsTableShortCol::operator()";
            std::cout << std::endl << "  Reference sum: " << sum_sht;
            std::cout << std::endl << "  Derived sum:   " << tot_sht << std::endl;
            throw;
        }
        std::cout << ".";

        // Check short table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i)
            tot_dbl += col_sht->real(i);
        if (!dequal(tot_dbl, double(sum_sht))) {
            std::cout << std::endl << "TEST ERROR: GFitsTableShortCol::real()";
            std::cout << std::endl << "  Reference sum: " << sum_sht;
            std::cout << std::endl << "  Derived sum:   " << tot_dbl << std::endl;
            throw;
        }
        std::cout << ".";

        // Check short table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i)
            tot_int += col_sht->integer(i);
        if (tot_int != sum_sht_int) {
            std::cout << std::endl << "TEST ERROR: GFitsTableShortCol::integer()";
            std::cout << std::endl << "  Reference sum: " << sum_sht_int;
            std::cout << std::endl << "  Derived sum:   " << tot_int << std::endl;
            throw;
        }
        std::cout << ".";

        // Check single precision table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i)
            tot_str += ":"+col_sht->string(i);
        if (tot_str != sum_sht_str) {
            std::cout << std::endl << "TEST ERROR: GFitsTableShortCol::string()";
            std::cout << std::endl << "  Reference sum: " << sum_sht_str;
            std::cout << std::endl << "  Derived sum:   " << tot_str << std::endl;
            throw;
        }
        std::cout << ".";

        //
        // ===== S H O R T  1 0 =====
        //

        // Get column
        GFitsTableShortCol* col_sht10 =
            (GFitsTableShortCol*)((GFitsTable*)fits.hdu(1))->column("SHORT10");

        // Check short table (operator access)
        tot_sht = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_sht += (*col_sht10)(i,j);
        }
        if (tot_sht != sum_sht10) {
            std::cout << std::endl << "TEST ERROR: GFitsTableShortCol::operator() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_sht10;
            std::cout << std::endl << "  Derived sum:   " << tot_sht << std::endl;
            throw;
        }
        std::cout << ".";

        // Check short table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_dbl += col_sht10->real(i,j);
        }
        if (!dequal(tot_dbl, sum_sht10)) {
            std::cout << std::endl << "TEST ERROR: GFitsTableShortCol::real() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_sht10;
            std::cout << std::endl << "  Derived sum:   " << tot_dbl << std::endl;
            throw;
        }
        std::cout << ".";

        // Check short table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_int += col_sht10->integer(i,j);
        }
        if (tot_int != sum_sht10_int) {
            std::cout << std::endl << "TEST ERROR: GFitsTableShortCol::integer() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_sht10_int;
            std::cout << std::endl << "  Derived sum:   " << tot_int << std::endl;
            throw;
        }
        std::cout << ".";

        // Check single precision table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_str += ":"+col_sht10->string(i,j);
        }
        if (tot_str != sum_sht10_str) {
            std::cout << std::endl << "TEST ERROR: GFitsTableShortCol::string() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_sht10_str;
            std::cout << std::endl << "  Derived sum:   " << tot_str << std::endl;
            throw;
        }
        std::cout << ".";
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to read GFitsTableShortCol tables."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ". ok." << std::endl;

}


/***************************************************************************
 * @brief Test unsigned short FITS binary table
 ***************************************************************************/
void test_bintable_ushort(void)
{
    // Set filename
    std::string filename = "test_bintable_ushort.fits";

    // Dump header
    std::cout << "Test GFitsTableUShortCol: ";

    // Remove FITS file
    std::string cmd = "rm -rf "+ filename;
    system(cmd.c_str());

    // Allocate reference sums
    unsigned short sum_sht       = 0;
    unsigned short sum_sht10     = 0;
    unsigned short sum_sht_int   = 0;
    unsigned short sum_sht10_int = 0;
    std::string    sum_sht_str;
    std::string    sum_sht10_str;

    // Build tables
    try {
        // Open FITS file
        GFits fits;
        fits.open(filename);

        // Set number of rows
        int nrows = 10;
        int nvec  = 10;

        // Initial control sums
        float          tot_flt = 0.0;
        double         tot_dbl = 0.0;
        int            tot_int = 0;
        long           tot_lng = 0;
        unsigned short tot_sht = 0;
        std::string    tot_str;

        //
        // ===== U S H O R T =====
        //

        // Set table
        GFitsTableUShortCol col_sht = GFitsTableUShortCol("USHORT", nrows);
        for (int i = 0; i < nrows; ++i) {
            unsigned short val_sht = (unsigned short)abs(1000.0 * cos(0.1*float(i)));
            col_sht(i)    = val_sht;
            sum_sht      += val_sht;
            sum_sht_int  += int(val_sht);
            sum_sht_str += ":"+str(val_sht);
        }
        std::cout << ".";

        // Check table (operator access)
        tot_sht = 0;
        for (int i = 0; i < nrows; ++i)
            tot_sht += col_sht(i);
        if (tot_sht != sum_sht) {
            std::cout << std::endl << "TEST ERROR: GFitsTableUShortCol::operator()";
            std::cout << std::endl << "  Reference sum: " << sum_sht;
            std::cout << std::endl << "  Derived sum:   " << tot_sht << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i)
            tot_dbl += col_sht.real(i);
        if (!dequal(tot_dbl, double(sum_sht))) {
            std::cout << std::endl << "TEST ERROR: GFitsTableUShortCol::real()";
            std::cout << std::endl << "  Reference sum: " << sum_sht;
            std::cout << std::endl << "  Derived sum:   " << tot_dbl << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i)
            tot_int += col_sht.integer(i);
        if (tot_int != sum_sht_int) {
            std::cout << std::endl << "TEST ERROR: GFitsTableUShortCol::integer()";
            std::cout << std::endl << "  Reference sum: " << sum_sht_int;
            std::cout << std::endl << "  Derived sum:   " << tot_int << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i)
            tot_str += ":"+col_sht.string(i);
        if (tot_str != sum_sht_str) {
            std::cout << std::endl << "TEST ERROR: GFitsTableUShortCol::string()";
            std::cout << std::endl << "  Reference sum: " << sum_sht_str;
            std::cout << std::endl << "  Derived sum:   " << tot_str << std::endl;
            throw;
        }
        std::cout << ".";

        //
        // ===== S H O R T  1 0 =====
        //

        // Set table
        GFitsTableUShortCol col_sht10 = GFitsTableUShortCol("USHORT10", nrows, nvec);
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j) {
                unsigned short val_sht  =
                  (unsigned short)(abs(100.0*cos(0.1*float(i))*cos(0.33*float(j))));
                col_sht10(i,j) = val_sht;
                sum_sht10     += val_sht;
                sum_sht10_int += int(val_sht);
                sum_sht10_str += ":"+str(val_sht);
            }
        }
        std::cout << ".";

        // Check table (operator access)
        tot_sht = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_sht += col_sht10(i,j);
        }
        if (tot_sht != sum_sht10) {
            std::cout << std::endl << "TEST ERROR: GFitsTableUShortCol::operator() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_sht10;
            std::cout << std::endl << "  Derived sum:   " << tot_sht << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_dbl += col_sht10.real(i,j);
        }
        if (!dequal(tot_dbl, sum_sht10)) {
            std::cout << std::endl << "TEST ERROR: GFitsTableUShortCol::real() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_sht10;
            std::cout << std::endl << "  Derived sum:   " << tot_dbl << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_int += col_sht10.integer(i,j);
        }
        if (tot_int != sum_sht10_int) {
            std::cout << std::endl << "TEST ERROR: GFitsTableUShortCol::integer() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_sht10_int;
            std::cout << std::endl << "  Derived sum:   " << tot_int << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_str += ":"+col_sht10.string(i,j);
        }
        if (tot_str != sum_sht10_str) {
            std::cout << std::endl << "TEST ERROR: GFitsTableUShortCol::string() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_sht10_str;
            std::cout << std::endl << "  Derived sum:   " << tot_str << std::endl;
            throw;
        }
        std::cout << ".";

        //
        // ===== W R I T E   T A B L E =====
        //

        // Build binary table
        GFitsBinTable table = GFitsBinTable(nrows);
        table.append_column(col_sht);
        table.insert_column(0, col_sht10);

        // Append to file
        fits.append(&table);

        // Save FITS file
        fits.save();
        std::cout << ".";

    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to build GFitsTableUShortCol tables."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";


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
        float          tot_flt = 0.0;
        unsigned short tot_sht = 0;
        double         tot_dbl = 0.0;
        int            tot_int = 0;
        std::string    tot_str;

        //
        // ===== U S H O R T =====
        //

        // Get column
        GFitsTableUShortCol* col_sht =
            (GFitsTableUShortCol*)((GFitsTable*)fits.hdu(1))->column("USHORT");

        // Check table (operator access)
        tot_sht = 0;
        for (int i = 0; i < nrows; ++i)
            tot_sht += (*col_sht)(i);
        if (tot_sht != sum_sht) {
            std::cout << std::endl << "TEST ERROR: GFitsTableUShortCol::operator()";
            std::cout << std::endl << "  Reference sum: " << sum_sht;
            std::cout << std::endl << "  Derived sum:   " << tot_sht << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i)
            tot_dbl += col_sht->real(i);
        if (!dequal(tot_dbl, double(sum_sht))) {
            std::cout << std::endl << "TEST ERROR: GFitsTableUShortCol::real()";
            std::cout << std::endl << "  Reference sum: " << sum_sht;
            std::cout << std::endl << "  Derived sum:   " << tot_dbl << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i)
            tot_int += col_sht->integer(i);
        if (tot_int != sum_sht_int) {
            std::cout << std::endl << "TEST ERROR: GFitsTableUShortCol::integer()";
            std::cout << std::endl << "  Reference sum: " << sum_sht_int;
            std::cout << std::endl << "  Derived sum:   " << tot_int << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i)
            tot_str += ":"+col_sht->string(i);
        if (tot_str != sum_sht_str) {
            std::cout << std::endl << "TEST ERROR: GFitsTableUShortCol::string()";
            std::cout << std::endl << "  Reference sum: " << sum_sht_str;
            std::cout << std::endl << "  Derived sum:   " << tot_str << std::endl;
            throw;
        }
        std::cout << ".";

        //
        // ===== U S H O R T  1 0 =====
        //

        // Get column
        GFitsTableUShortCol* col_sht10 =
            (GFitsTableUShortCol*)((GFitsTable*)fits.hdu(1))->column("USHORT10");

        // Check table (operator access)
        tot_sht = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_sht += (*col_sht10)(i,j);
        }
        if (tot_sht != sum_sht10) {
            std::cout << std::endl << "TEST ERROR: GFitsTableUShortCol::operator() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_sht10;
            std::cout << std::endl << "  Derived sum:   " << tot_sht << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_dbl += col_sht10->real(i,j);
        }
        if (!dequal(tot_dbl, sum_sht10)) {
            std::cout << std::endl << "TEST ERROR: GFitsTableUShortCol::real() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_sht10;
            std::cout << std::endl << "  Derived sum:   " << tot_dbl << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_int += col_sht10->integer(i,j);
        }
        if (tot_int != sum_sht10_int) {
            std::cout << std::endl << "TEST ERROR: GFitsTableUShortCol::integer() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_sht10_int;
            std::cout << std::endl << "  Derived sum:   " << tot_int << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_str += ":"+col_sht10->string(i,j);
        }
        if (tot_str != sum_sht10_str) {
            std::cout << std::endl << "TEST ERROR: GFitsTableUShortCol::string() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_sht10_str;
            std::cout << std::endl << "  Derived sum:   " << tot_str << std::endl;
            throw;
        }
        std::cout << ".";
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to read GFitsTableUShortCol tables."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ". ok." << std::endl;

}


/***************************************************************************
 * @brief Test long FITS binary table
 ***************************************************************************/
void test_bintable_long(void)
{
    // Set filename
    std::string filename = "test_bintable_long.fits";

    // Dump header
    std::cout << "Test GFitsTableLongCol: ";
    
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
        GFitsTableLongCol col_lng = GFitsTableLongCol("LONG", nrows);
        for (int i = 0; i < nrows; ++i) {
            long val_lng  = long(100000.0 * cos(0.1*float(i)));
            col_lng(i)    = val_lng;
            sum_lng      += val_lng;
            sum_lng_int  += int(val_lng);
            sum_lng_str += ":"+str((int)val_lng);
        }
        std::cout << ".";

        // Check long table (operator access)
        tot_lng = 0;
        for (int i = 0; i < nrows; ++i)
            tot_lng += col_lng(i);
        if (tot_lng != sum_lng) {
            std::cout << std::endl << "TEST ERROR: GFitsTableLongCol::operator()";
            std::cout << std::endl << "  Reference sum: " << sum_lng;
            std::cout << std::endl << "  Derived sum:   " << tot_lng << std::endl;
            throw;
        }
        std::cout << ".";

        // Check long table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i)
            tot_dbl += col_lng.real(i);
        if (!dequal(tot_dbl, double(sum_lng))) {
            std::cout << std::endl << "TEST ERROR: GFitsTableLongCol::real()";
            std::cout << std::endl << "  Reference sum: " << sum_lng;
            std::cout << std::endl << "  Derived sum:   " << tot_dbl << std::endl;
            throw;
        }
        std::cout << ".";

        // Check long table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i)
            tot_int += col_lng.integer(i);
        if (tot_int != sum_lng_int) {
            std::cout << std::endl << "TEST ERROR: GFitsTableLongCol::integer()";
            std::cout << std::endl << "  Reference sum: " << sum_lng_int;
            std::cout << std::endl << "  Derived sum:   " << tot_int << std::endl;
            throw;
        }
        std::cout << ".";

        // Check long table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i)
            tot_str += ":"+col_lng.string(i);
        if (tot_str != sum_lng_str) {
            std::cout << std::endl << "TEST ERROR: GFitsTableLongCol::string()";
            std::cout << std::endl << "  Reference sum: " << sum_lng_str;
            std::cout << std::endl << "  Derived sum:   " << tot_str << std::endl;
            throw;
        }
        std::cout << ".";

        //
        // ===== L O N G  1 0 =====
        //

        // Set long table
        GFitsTableLongCol col_lng10 = GFitsTableLongCol("LONG10", nrows, nvec);
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j) {
                long val_lng   = long(1000.0*cos(0.1*float(i))*
                                             cos(0.33*float(j)));
                col_lng10(i,j) = val_lng;
                sum_lng10     += val_lng;
                sum_lng10_int += int(val_lng);
                sum_lng10_str += ":"+str((int)val_lng);
            }
        }
        std::cout << ".";

        // Check long table (operator access)
        tot_lng = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_lng += col_lng10(i,j);
        }
        if (tot_lng != sum_lng10) {
            std::cout << std::endl << "TEST ERROR: GFitsTableLongCol::operator() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_lng10;
            std::cout << std::endl << "  Derived sum:   " << tot_lng << std::endl;
            throw;
        }
        std::cout << ".";

        // Check long table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_dbl += col_lng10.real(i,j);
        }
        if (!dequal(tot_dbl, sum_lng10)) {
            std::cout << std::endl << "TEST ERROR: GFitsTableLongCol::real() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_lng10;
            std::cout << std::endl << "  Derived sum:   " << tot_dbl << std::endl;
            throw;
        }
        std::cout << ".";

        // Check long table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_int += col_lng10.integer(i,j);
        }
        if (tot_int != sum_lng10_int) {
            std::cout << std::endl << "TEST ERROR: GFitsTableLongCol::integer() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_lng10_int;
            std::cout << std::endl << "  Derived sum:   " << tot_int << std::endl;
            throw;
        }
        std::cout << ".";

        // Check single precision table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_str += ":"+col_lng10.string(i,j);
        }
        if (tot_str != sum_lng10_str) {
            std::cout << std::endl << "TEST ERROR: GFitsTableLongCol::string() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_lng10_str;
            std::cout << std::endl << "  Derived sum:   " << tot_str << std::endl;
            throw;
        }
        std::cout << ".";

        //
        // ===== W R I T E   T A B L E =====
        //

        // Build binary table
        GFitsBinTable table = GFitsBinTable(nrows);
        table.append_column(col_lng);
        table.insert_column(99, col_lng10);

        // Append to file
        fits.append(&table);

        // Save FITS file
        fits.save();
        std::cout << ".";
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to build GFitsTableLongCol tables."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";


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
        GFitsTableLongCol* col_lng =
            (GFitsTableLongCol*)((GFitsTable*)fits.hdu(1))->column("LONG");

        // Check long table (operator access)
        tot_lng = 0;
        for (int i = 0; i < nrows; ++i)
            tot_lng += (*col_lng)(i);
        if (tot_lng != sum_lng) {
            std::cout << std::endl << "TEST ERROR: GFitsTableLongCol::operator()";
            std::cout << std::endl << "  Reference sum: " << sum_lng;
            std::cout << std::endl << "  Derived sum:   " << tot_lng << std::endl;
            throw;
        }
        std::cout << ".";

        // Check long table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i)
            tot_dbl += col_lng->real(i);
        if (!dequal(tot_dbl, double(sum_lng))) {
            std::cout << std::endl << "TEST ERROR: GFitsTableLongCol::real()";
            std::cout << std::endl << "  Reference sum: " << sum_lng;
            std::cout << std::endl << "  Derived sum:   " << tot_dbl << std::endl;
            throw;
        }
        std::cout << ".";

        // Check long table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i)
            tot_int += col_lng->integer(i);
        if (tot_int != sum_lng_int) {
            std::cout << std::endl << "TEST ERROR: GFitsTableLongCol::integer()";
            std::cout << std::endl << "  Reference sum: " << sum_lng_int;
            std::cout << std::endl << "  Derived sum:   " << tot_int << std::endl;
            throw;
        }
        std::cout << ".";

        // Check long table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i)
            tot_str += ":"+col_lng->string(i);
        if (tot_str != sum_lng_str) {
            std::cout << std::endl << "TEST ERROR: GFitsTableLongCol::string()";
            std::cout << std::endl << "  Reference sum: " << sum_lng_str;
            std::cout << std::endl << "  Derived sum:   " << tot_str << std::endl;
            throw;
        }
        std::cout << ".";

        //
        // ===== L O N G  1 0 =====
        //

        // Get column
        GFitsTableLongCol* col_lng10 =
            (GFitsTableLongCol*)((GFitsTable*)fits.hdu(1))->column("LONG10");

        // Check long table (operator access)
        tot_lng = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_lng += (*col_lng10)(i,j);
        }
        if (tot_lng != sum_lng10) {
            std::cout << std::endl << "TEST ERROR: GFitsTableLongCol::operator() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_lng10;
            std::cout << std::endl << "  Derived sum:   " << tot_lng << std::endl;
            throw;
        }
        std::cout << ".";

        // Check long table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_dbl += col_lng10->real(i,j);
        }
        if (!dequal(tot_dbl, sum_lng10)) {
            std::cout << std::endl << "TEST ERROR: GFitsTableLongCol::real() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_lng10;
            std::cout << std::endl << "  Derived sum:   " << tot_dbl << std::endl;
            throw;
        }
        std::cout << ".";

        // Check long table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_int += col_lng10->integer(i,j);
        }
        if (tot_int != sum_lng10_int) {
            std::cout << std::endl << "TEST ERROR: GFitsTableLongCol::integer() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_lng10_int;
            std::cout << std::endl << "  Derived sum:   " << tot_int << std::endl;
            throw;
        }
        std::cout << ".";

        // Check single precision table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_str += ":"+col_lng10->string(i,j);
        }
        if (tot_str != sum_lng10_str) {
            std::cout << std::endl << "TEST ERROR: GFitsTableLongCol::string() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_lng10_str;
            std::cout << std::endl << "  Derived sum:   " << tot_str << std::endl;
            throw;
        }
        std::cout << ".";
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to read GFitsTableLongCol tables."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ". ok." << std::endl;

}


/***************************************************************************
 * @brief Test long long FITS binary table
 ***************************************************************************/
void test_bintable_longlong(void)
{
    // Set filename
    std::string filename = "test_bintable_longlong.fits";

    // Dump header
    std::cout << "Test GFitsTableLongLongCol: ";
    
    // Remove FITS file
    std::string cmd = "rm -rf "+ filename;
    system(cmd.c_str());

    // Allocate reference sums
    long long   sum_lng       = 0;
    long long   sum_lng10     = 0;
    long long   sum_lng_int   = 0;
    long long   sum_lng10_int = 0;
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
        long long   tot_lng = 0;
        short       tot_sht = 0;
        std::string tot_str;

        //
        // ===== L O N G  L O N G=====
        //

        // Set table
        GFitsTableLongLongCol col_lng = GFitsTableLongLongCol("LONGLONG", nrows);
        for (int i = 0; i < nrows; ++i) {
            long val_lng  = long(100000.0 * cos(0.1*float(i)));
            col_lng(i)    = val_lng;
            sum_lng      += val_lng;
            sum_lng_int  += int(val_lng);
            sum_lng_str += ":"+str((int)val_lng);
        }
        std::cout << ".";

        // Check table (operator access)
        tot_lng = 0;
        for (int i = 0; i < nrows; ++i)
            tot_lng += col_lng(i);
        if (tot_lng != sum_lng) {
            std::cout << std::endl << "TEST ERROR: GFitsTableLongLongCol::operator()";
            std::cout << std::endl << "  Reference sum: " << sum_lng;
            std::cout << std::endl << "  Derived sum:   " << tot_lng << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i)
            tot_dbl += col_lng.real(i);
        if (!dequal(tot_dbl, double(sum_lng))) {
            std::cout << std::endl << "TEST ERROR: GFitsTableLongLongCol::real()";
            std::cout << std::endl << "  Reference sum: " << sum_lng;
            std::cout << std::endl << "  Derived sum:   " << tot_dbl << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i)
            tot_int += col_lng.integer(i);
        if (tot_int != sum_lng_int) {
            std::cout << std::endl << "TEST ERROR: GFitsTableLongLongCol::integer()";
            std::cout << std::endl << "  Reference sum: " << sum_lng_int;
            std::cout << std::endl << "  Derived sum:   " << tot_int << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i)
            tot_str += ":"+col_lng.string(i);
        if (tot_str != sum_lng_str) {
            std::cout << std::endl << "TEST ERROR: GFitsTableLongLongCol::string()";
            std::cout << std::endl << "  Reference sum: " << sum_lng_str;
            std::cout << std::endl << "  Derived sum:   " << tot_str << std::endl;
            throw;
        }
        std::cout << ".";

        //
        // ===== L O N G  L O N G 1 0 =====
        //

        // Set long table
        GFitsTableLongLongCol col_lng10 = GFitsTableLongLongCol("LONGLONG10", nrows, nvec);
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j) {
                long val_lng   = long(1000.0*cos(0.1*float(i))*
                                             cos(0.33*float(j)));
                col_lng10(i,j) = val_lng;
                sum_lng10     += val_lng;
                sum_lng10_int += int(val_lng);
                sum_lng10_str += ":"+str((int)val_lng);
            }
        }
        std::cout << ".";

        // Check table (operator access)
        tot_lng = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_lng += col_lng10(i,j);
        }
        if (tot_lng != sum_lng10) {
            std::cout << std::endl << "TEST ERROR: GFitsTableLongLongCol::operator() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_lng10;
            std::cout << std::endl << "  Derived sum:   " << tot_lng << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_dbl += col_lng10.real(i,j);
        }
        if (!dequal(tot_dbl, sum_lng10)) {
            std::cout << std::endl << "TEST ERROR: GFitsTableLongLongCol::real() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_lng10;
            std::cout << std::endl << "  Derived sum:   " << tot_dbl << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_int += col_lng10.integer(i,j);
        }
        if (tot_int != sum_lng10_int) {
            std::cout << std::endl << "TEST ERROR: GFitsTableLongLongCol::integer() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_lng10_int;
            std::cout << std::endl << "  Derived sum:   " << tot_int << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_str += ":"+col_lng10.string(i,j);
        }
        if (tot_str != sum_lng10_str) {
            std::cout << std::endl << "TEST ERROR: GFitsTableLongLongCol::string() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_lng10_str;
            std::cout << std::endl << "  Derived sum:   " << tot_str << std::endl;
            throw;
        }
        std::cout << ".";

        //
        // ===== W R I T E   T A B L E =====
        //

        // Build binary table
        GFitsBinTable table = GFitsBinTable(nrows);
        table.append_column(col_lng);
        table.insert_column(99, col_lng10);

        // Append to file
        fits.append(&table);

        // Save FITS file
        fits.save();
        std::cout << ".";
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to build GFitsTableLongLongCol tables."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";


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
        long long   tot_lng = 0;
        std::string tot_str;

        //
        // ===== L O N G L O N G =====
        //

        // Get column
        GFitsTableLongLongCol* col_lng =
            (GFitsTableLongLongCol*)((GFitsTable*)fits.hdu(1))->column("LONGLONG");

        // Check table (operator access)
        tot_lng = 0;
        for (int i = 0; i < nrows; ++i)
            tot_lng += (*col_lng)(i);
        if (tot_lng != sum_lng) {
            std::cout << std::endl << "TEST ERROR: GFitsTableLongLongCol::operator()";
            std::cout << std::endl << "  Reference sum: " << sum_lng;
            std::cout << std::endl << "  Derived sum:   " << tot_lng << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i)
            tot_dbl += col_lng->real(i);
        if (!dequal(tot_dbl, double(sum_lng))) {
            std::cout << std::endl << "TEST ERROR: GFitsTableLongLongCol::real()";
            std::cout << std::endl << "  Reference sum: " << sum_lng;
            std::cout << std::endl << "  Derived sum:   " << tot_dbl << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i)
            tot_int += col_lng->integer(i);
        if (tot_int != sum_lng_int) {
            std::cout << std::endl << "TEST ERROR: GFitsTableLongLongCol::integer()";
            std::cout << std::endl << "  Reference sum: " << sum_lng_int;
            std::cout << std::endl << "  Derived sum:   " << tot_int << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i)
            tot_str += ":"+col_lng->string(i);
        if (tot_str != sum_lng_str) {
            std::cout << std::endl << "TEST ERROR: GFitsTableLongLongCol::string()";
            std::cout << std::endl << "  Reference sum: " << sum_lng_str;
            std::cout << std::endl << "  Derived sum:   " << tot_str << std::endl;
            throw;
        }
        std::cout << ".";

        //
        // ===== L O N G  L O N G 1 0 =====
        //

        // Get column
        GFitsTableLongLongCol* col_lng10 =
            (GFitsTableLongLongCol*)((GFitsTable*)fits.hdu(1))->column("LONGLONG10");

        // Check table (operator access)
        tot_lng = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_lng += (*col_lng10)(i,j);
        }
        if (tot_lng != sum_lng10) {
            std::cout << std::endl << "TEST ERROR: GFitsTableLongLongCol::operator() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_lng10;
            std::cout << std::endl << "  Derived sum:   " << tot_lng << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_dbl += col_lng10->real(i,j);
        }
        if (!dequal(tot_dbl, sum_lng10)) {
            std::cout << std::endl << "TEST ERROR: GFitsTableLongLongCol::real() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_lng10;
            std::cout << std::endl << "  Derived sum:   " << tot_dbl << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_int += col_lng10->integer(i,j);
        }
        if (tot_int != sum_lng10_int) {
            std::cout << std::endl << "TEST ERROR: GFitsTableLongLongCol::integer() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_lng10_int;
            std::cout << std::endl << "  Derived sum:   " << tot_int << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_str += ":"+col_lng10->string(i,j);
        }
        if (tot_str != sum_lng10_str) {
            std::cout << std::endl << "TEST ERROR: GFitsTableLongLongCol::string() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_lng10_str;
            std::cout << std::endl << "  Derived sum:   " << tot_str << std::endl;
            throw;
        }
        std::cout << ".";
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to read GFitsTableLongLongCol tables."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ". ok." << std::endl;

}


/***************************************************************************
 * @brief Test unsigned long FITS binary table
 ***************************************************************************/
void test_bintable_ulong(void)
{
    // Set filename
    std::string filename = "test_bintable_ulong.fits";

    // Dump header
    std::cout << "Test GFitsTableULongCol: ";

    // Remove FITS file
    std::string cmd = "rm -rf "+ filename;
    system(cmd.c_str());

    // Allocate reference sums
    unsigned long sum_lng       = 0;
    unsigned long sum_lng10     = 0;
    unsigned long sum_lng_int   = 0;
    unsigned long sum_lng10_int = 0;
    std::string   sum_lng_str;
    std::string   sum_lng10_str;

    // Build tables
    try {
        // Open FITS file
        GFits fits;
        fits.open(filename);

        // Set number of rows
        int nrows = 10;
        int nvec  = 10;

        // Initial control sums
        float         tot_flt = 0.0;
        double        tot_dbl = 0.0;
        int           tot_int = 0;
        unsigned long tot_lng = 0;
        short         tot_sht = 0;
        std::string   tot_str;

        //
        // ===== U L O N G=====
        //

        // Set table
        GFitsTableULongCol col_lng = GFitsTableULongCol("ULONG", nrows);
        for (int i = 0; i < nrows; ++i) {
            unsigned long val_lng  = (unsigned long)abs(100000.0 * cos(0.1*float(i)));
            col_lng(i)    = val_lng;
            sum_lng      += val_lng;
            sum_lng_int  += int(val_lng);
            sum_lng_str += ":"+str((int)val_lng);
        }
        std::cout << ".";

        // Check table (operator access)
        tot_lng = 0;
        for (int i = 0; i < nrows; ++i)
            tot_lng += col_lng(i);
        if (tot_lng != sum_lng) {
            std::cout << std::endl << "TEST ERROR: GFitsTableULongCol::operator()";
            std::cout << std::endl << "  Reference sum: " << sum_lng;
            std::cout << std::endl << "  Derived sum:   " << tot_lng << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i)
            tot_dbl += col_lng.real(i);
        if (!dequal(tot_dbl, double(sum_lng))) {
            std::cout << std::endl << "TEST ERROR: GFitsTableULongCol::real()";
            std::cout << std::endl << "  Reference sum: " << sum_lng;
            std::cout << std::endl << "  Derived sum:   " << tot_dbl << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i)
            tot_int += col_lng.integer(i);
        if (tot_int != sum_lng_int) {
            std::cout << std::endl << "TEST ERROR: GFitsTableULongCol::integer()";
            std::cout << std::endl << "  Reference sum: " << sum_lng_int;
            std::cout << std::endl << "  Derived sum:   " << tot_int << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i)
            tot_str += ":"+col_lng.string(i);
        if (tot_str != sum_lng_str) {
            std::cout << std::endl << "TEST ERROR: GFitsTableULongCol::string()";
            std::cout << std::endl << "  Reference sum: " << sum_lng_str;
            std::cout << std::endl << "  Derived sum:   " << tot_str << std::endl;
            throw;
        }
        std::cout << ".";

        //
        // ===== U L O N G 1 0 =====
        //

        // Set long table
        GFitsTableULongCol col_lng10 = GFitsTableULongCol("ULONG10", nrows, nvec);
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j) {
                unsigned long val_lng = (unsigned long)abs(1000.0*cos(0.1*float(i))*
                                                           cos(0.33*float(j)));
                col_lng10(i,j) = val_lng;
                sum_lng10     += val_lng;
                sum_lng10_int += int(val_lng);
                sum_lng10_str += ":"+str((int)val_lng);
            }
        }
        std::cout << ".";

        // Check table (operator access)
        tot_lng = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_lng += col_lng10(i,j);
        }
        if (tot_lng != sum_lng10) {
            std::cout << std::endl << "TEST ERROR: GFitsTableULongCol::operator() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_lng10;
            std::cout << std::endl << "  Derived sum:   " << tot_lng << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_dbl += col_lng10.real(i,j);
        }
        if (!dequal(tot_dbl, sum_lng10)) {
            std::cout << std::endl << "TEST ERROR: GFitsTableULongCol::real() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_lng10;
            std::cout << std::endl << "  Derived sum:   " << tot_dbl << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_int += col_lng10.integer(i,j);
        }
        if (tot_int != sum_lng10_int) {
            std::cout << std::endl << "TEST ERROR: GFitsTableULongCol::integer() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_lng10_int;
            std::cout << std::endl << "  Derived sum:   " << tot_int << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_str += ":"+col_lng10.string(i,j);
        }
        if (tot_str != sum_lng10_str) {
            std::cout << std::endl << "TEST ERROR: GFitsTableULongCol::string() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_lng10_str;
            std::cout << std::endl << "  Derived sum:   " << tot_str << std::endl;
            throw;
        }
        std::cout << ".";

        //
        // ===== W R I T E   T A B L E =====
        //

        // Build binary table
        GFitsBinTable table = GFitsBinTable(nrows);
        table.append_column(col_lng);
        table.insert_column(99, col_lng10);

        // Append to file
        fits.append(&table);

        // Save FITS file
        fits.save();
        std::cout << ".";
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to build GFitsTableULongCol tables."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";


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
        float         tot_flt = 0.0;
        double        tot_dbl = 0.0;
        int           tot_int = 0;
        unsigned long tot_lng = 0;
        std::string   tot_str;

        //
        // ===== U L O N G =====
        //

        // Get column
        GFitsTableULongCol* col_lng =
            (GFitsTableULongCol*)((GFitsTable*)fits.hdu(1))->column("ULONG");

        // Check table (operator access)
        tot_lng = 0;
        for (int i = 0; i < nrows; ++i)
            tot_lng += (*col_lng)(i);
        if (tot_lng != sum_lng) {
            std::cout << std::endl << "TEST ERROR: GFitsTableULongCol::operator()";
            std::cout << std::endl << "  Reference sum: " << sum_lng;
            std::cout << std::endl << "  Derived sum:   " << tot_lng << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i)
            tot_dbl += col_lng->real(i);
        if (!dequal(tot_dbl, double(sum_lng))) {
            std::cout << std::endl << "TEST ERROR: GFitsTableULongCol::real()";
            std::cout << std::endl << "  Reference sum: " << sum_lng;
            std::cout << std::endl << "  Derived sum:   " << tot_dbl << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i)
            tot_int += col_lng->integer(i);
        if (tot_int != sum_lng_int) {
            std::cout << std::endl << "TEST ERROR: GFitsTableULongCol::integer()";
            std::cout << std::endl << "  Reference sum: " << sum_lng_int;
            std::cout << std::endl << "  Derived sum:   " << tot_int << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i)
            tot_str += ":"+col_lng->string(i);
        if (tot_str != sum_lng_str) {
            std::cout << std::endl << "TEST ERROR: GFitsTableULongCol::string()";
            std::cout << std::endl << "  Reference sum: " << sum_lng_str;
            std::cout << std::endl << "  Derived sum:   " << tot_str << std::endl;
            throw;
        }
        std::cout << ".";

        //
        // ===== U L O N G 1 0 =====
        //

        // Get column
        GFitsTableULongCol* col_lng10 =
            (GFitsTableULongCol*)((GFitsTable*)fits.hdu(1))->column("ULONG10");

        // Check table (operator access)
        tot_lng = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_lng += (*col_lng10)(i,j);
        }
        if (tot_lng != sum_lng10) {
            std::cout << std::endl << "TEST ERROR: GFitsTableULongCol::operator() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_lng10;
            std::cout << std::endl << "  Derived sum:   " << tot_lng << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_dbl += col_lng10->real(i,j);
        }
        if (!dequal(tot_dbl, sum_lng10)) {
            std::cout << std::endl << "TEST ERROR: GFitsTableULongCol::real() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_lng10;
            std::cout << std::endl << "  Derived sum:   " << tot_dbl << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_int += col_lng10->integer(i,j);
        }
        if (tot_int != sum_lng10_int) {
            std::cout << std::endl << "TEST ERROR: GFitsTableULongCol::integer() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_lng10_int;
            std::cout << std::endl << "  Derived sum:   " << tot_int << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_str += ":"+col_lng10->string(i,j);
        }
        if (tot_str != sum_lng10_str) {
            std::cout << std::endl << "TEST ERROR: GFitsTableULongCol::string() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_lng10_str;
            std::cout << std::endl << "  Derived sum:   " << tot_str << std::endl;
            throw;
        }
        std::cout << ".";
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to read GFitsTableULongCol tables."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ". ok." << std::endl;

}


/***************************************************************************
 * @brief Test string FITS binary table
 ***************************************************************************/
void test_bintable_string(void)
{
    // Set filename
    std::string filename = "test_bintable_string.fits";

    // Dump header
    std::cout << "Test GFitsTableStringCol: ";

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
        GFitsTableStringCol col_str = GFitsTableStringCol("STRING", nrows, 20);
        for (int i = 0; i < nrows; ++i) {
            double val_dbl = 100.0 * cos(0.1*double(i));
            col_str(i)   = str(val_dbl);
            sum_str     += ":" + str(val_dbl);
            sum_str_dbl += val_dbl;
            sum_str_int += int(val_dbl);
        }
        std::cout << ".";

        // Check string table (operator access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i)
            tot_str += ":" + col_str(i);
        if (tot_str != sum_str) {
            std::cout << std::endl << "TEST ERROR: GFitsTableStringCol::operator()";
            std::cout << std::endl << "  Reference sum: " << sum_str;
            std::cout << std::endl << "  Derived sum:   " << tot_str << std::endl;
            throw;
        }
        std::cout << ".";

        // Check string table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i)
            tot_dbl += col_str.real(i);
        if (!equal(tot_dbl, sum_str_dbl, 1.0e-2)) {
            std::cout << std::endl << "TEST ERROR: GFitsTableStringCol::real()";
            std::cout << std::endl << "  Reference sum: " << sum_str_dbl;
            std::cout << std::endl << "  Derived sum:   " << tot_dbl << std::endl;
            throw;
        }
        std::cout << ".";

        // Check string table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i)
            tot_int += col_str.integer(i);
        if (tot_int != sum_str_int) {
            std::cout << std::endl << "TEST ERROR: GFitsTableStringCol::integer()";
            std::cout << std::endl << "  Reference sum: " << sum_str_int;
            std::cout << std::endl << "  Derived sum:   " << tot_int << std::endl;
            throw;
        }
        std::cout << ".";

        // Check string table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i)
            tot_str += ":"+col_str.string(i);
        if (tot_str != sum_str) {
            std::cout << std::endl << "TEST ERROR: GFitsTableStringCol::string()";
            std::cout << std::endl << "  Reference sum: " << sum_str;
            std::cout << std::endl << "  Derived sum:   " << tot_str << std::endl;
            throw;
        }
        std::cout << ".";

        //
        // ===== S T R I N G   1 0 =====
        //

        // Set string table
        GFitsTableStringCol col_str10 = GFitsTableStringCol("STRING10", nrows, 20, nvec);
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j) {
                double val_dbl = 100.0 * cos(0.1*double(i)) * cos(0.3*double(j));
                col_str10(i,j)  = str(val_dbl);
                sum_str10      += ":" + str(val_dbl);
                sum_str10_dbl  += val_dbl;
                sum_str10_int  += int(val_dbl);
            }
        }
        std::cout << ".";

        // Check string table (operator access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_str += ":" + col_str10(i,j);
        }
        if (tot_str != sum_str10) {
            std::cout << std::endl << "TEST ERROR: GFitsTableStringCol::operator() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_str10;
            std::cout << std::endl << "  Derived sum:   " << tot_str << std::endl;
            throw;
        }
        std::cout << ".";

        // Check string table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_dbl += col_str10.real(i,j);
        }
        if (!equal(tot_dbl, sum_str10_dbl, 1.0e-2)) {
            std::cout << std::endl << "TEST ERROR: GFitsTableStringCol::real() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_str10_dbl;
            std::cout << std::endl << "  Derived sum:   " << tot_dbl << std::endl;
            throw;
        }
        std::cout << ".";

        // Check string table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_int += col_str10.integer(i,j);
        }
        if (tot_int != sum_str10_int) {
            std::cout << std::endl << "TEST ERROR: GFitsTableStringCol::integer() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_str10_int;
            std::cout << std::endl << "  Derived sum:   " << tot_int << std::endl;
            throw;
        }
        std::cout << ".";

        // Check string table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_str += ":"+col_str10.string(i,j);
        }
        if (tot_str != sum_str10) {
            std::cout << std::endl << "TEST ERROR: GFitsTableStringCol::string() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_str10;
            std::cout << std::endl << "  Derived sum:   " << tot_str << std::endl;
            throw;
        }
        std::cout << ".";

        //
        // ===== W R I T E   T A B L E =====
        //

        // Build binary table
        GFitsBinTable table = GFitsBinTable(nrows);
        table.append_column(col_str);
        table.append_column(col_str10);

        // Append to file
        fits.append(&table);

        // Save FITS file
        fits.save();
        std::cout << ".";

    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to build GFitsTableStringCol tables."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";

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
        GFitsTableStringCol* col_str =
            (GFitsTableStringCol*)((GFitsTable*)fits.hdu(1))->column("STRING");

        // Check string table (operator access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i)
            tot_str += ":" + (*col_str)(i);
        if (tot_str != sum_str) {
            std::cout << std::endl << "TEST ERROR: GFitsTableStringCol::operator()";
            std::cout << std::endl << "  Reference sum: " << sum_str;
            std::cout << std::endl << "  Derived sum:   " << tot_str << std::endl;
            throw;
        }
        std::cout << ".";

        // Check string table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i)
            tot_dbl += col_str->real(i);
        if (!equal(tot_dbl, sum_str_dbl, 1.0e-2)) {
            std::cout << std::endl << "TEST ERROR: GFitsTableStringCol::real()";
            std::cout << std::endl << "  Reference sum: " << sum_str_dbl;
            std::cout << std::endl << "  Derived sum:   " << tot_dbl << std::endl;
            throw;
        }
        std::cout << ".";

        // Check string table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i)
            tot_int += col_str->integer(i);
        if (tot_int != sum_str_int) {
            std::cout << std::endl << "TEST ERROR: GFitsTableStringCol::integer()";
            std::cout << std::endl << "  Reference sum: " << sum_str_int;
            std::cout << std::endl << "  Derived sum:   " << tot_int << std::endl;
            throw;
        }
        std::cout << ".";

        // Check string table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i)
            tot_str += ":"+col_str->string(i);
        if (tot_str != sum_str) {
            std::cout << std::endl << "TEST ERROR: GFitsTableStringCol::string()";
            std::cout << std::endl << "  Reference sum: " << sum_str;
            std::cout << std::endl << "  Derived sum:   " << tot_str << std::endl;
            throw;
        }
        std::cout << ".";

        //
        // ===== S T R I N G   1 0 =====
        //

        // Get column
        GFitsTableStringCol* col_str10 =
            (GFitsTableStringCol*)((GFitsTable*)fits.hdu(1))->column("STRING10");

        // Check table (operator access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_str += ":" + (*col_str10)(i,j);
        }
        if (tot_str != sum_str10) {
            std::cout << std::endl << "TEST ERROR: GFitsTableStringCol::operator() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_str10;
            std::cout << std::endl << "  Derived sum:   " << tot_str << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_dbl += col_str10->real(i,j);
        }
        if (!equal(tot_dbl, sum_str10_dbl, 1.0e-2)) {
            std::cout << std::endl << "TEST ERROR: GFitsTableStringCol::real() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_str10_dbl;
            std::cout << std::endl << "  Derived sum:   " << tot_dbl << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_int += col_str10->integer(i,j);
        }
        if (tot_int != sum_str10_int) {
            std::cout << std::endl << "TEST ERROR: GFitsTableStringCol::integer() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_str10_int;
            std::cout << std::endl << "  Derived sum:   " << tot_int << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_str += ":"+col_str10->string(i,j);
        }
        if (tot_str != sum_str10) {
            std::cout << std::endl << "TEST ERROR: GFitsTableStringCol::string() - 10";
            std::cout << std::endl << "  Reference sum: " << sum_str10;
            std::cout << std::endl << "  Derived sum:   " << tot_str << std::endl;
            throw;
        }
        std::cout << ".";

    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to read GFitsTableStringCol tables."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ". ok." << std::endl;

}


/***************************************************************************
 * @brief Test boolean FITS binary table
 ***************************************************************************/
void test_bintable_logical(void)
{
    // Set filename
    std::string filename = "test_bintable_logical.fits";

    // Dump header
    std::cout << "Test GFitsTableBoolCol: ";

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
        GFitsTableBoolCol col = GFitsTableBoolCol("LOGICAL", nrows);
        for (int i = 0; i < nrows; ++i) {
            bool val  = (i % 2) ? true : false;
            col(i)    = val;
            sum      += val;
            sum_int  += int(val);
            if (val)
                sum_str += ":T";
            else
                sum_str += ":F";
        }
        std::cout << ".";

        // Check table (operator access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i)
            tot_int += col(i);
        if (tot_int != sum) {
            std::cout << std::endl << "TEST ERROR: GFitsTableBoolCol::operator()";
            std::cout << std::endl << "  Reference sum: " << sum;
            std::cout << std::endl << "  Derived sum:   " << tot_int << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i)
            tot_dbl += col.real(i);
        if (!dequal(tot_dbl, double(sum))) {
            std::cout << std::endl << "TEST ERROR: GFitsTableBoolCol::real()";
            std::cout << std::endl << "  Reference sum: " << sum;
            std::cout << std::endl << "  Derived sum:   " << tot_dbl << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i)
            tot_int += col.integer(i);
        if (tot_int != sum_int) {
            std::cout << std::endl << "TEST ERROR: GFitsTableBoolCol::integer()";
            std::cout << std::endl << "  Reference sum: " << sum_int;
            std::cout << std::endl << "  Derived sum:   " << tot_int << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i)
            tot_str += ":"+col.string(i);
        if (tot_str != sum_str) {
            std::cout << std::endl << "TEST ERROR: GFitsTableBoolCol::string()";
            std::cout << std::endl << "  Reference sum: " << sum_str;
            std::cout << std::endl << "  Derived sum:   " << tot_str << std::endl;
            throw;
        }
        std::cout << ".";

        //
        // ===== L O G I C A L  1 0 =====
        //

        // Set table
        GFitsTableBoolCol col10 = GFitsTableBoolCol("LOGICAL10", nrows, nvec);
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
        std::cout << ".";

        // Check table (operator access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_int += col10(i,j);
        }
        if (tot_int != sum10) {
            std::cout << std::endl << "TEST ERROR: GFitsTableBoolCol::operator() - 10";
            std::cout << std::endl << "  Reference sum: " << sum10;
            std::cout << std::endl << "  Derived sum:   " << tot_int << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_dbl += col10.real(i,j);
        }
        if (!dequal(tot_dbl, sum10)) {
            std::cout << std::endl << "TEST ERROR: GFitsTableBoolCol::real() - 10";
            std::cout << std::endl << "  Reference sum: " << sum10;
            std::cout << std::endl << "  Derived sum:   " << tot_dbl << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_int += col10.integer(i,j);
        }
        if (tot_int != sum10_int) {
            std::cout << std::endl << "TEST ERROR: GFitsTableBoolCol::integer() - 10";
            std::cout << std::endl << "  Reference sum: " << sum10_int;
            std::cout << std::endl << "  Derived sum:   " << tot_int << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_str += ":"+col10.string(i,j);
        }
        if (tot_str != sum10_str) {
            std::cout << std::endl << "TEST ERROR: GFitsTableBoolCol::string() - 10";
            std::cout << std::endl << "  Reference sum: " << sum10_str;
            std::cout << std::endl << "  Derived sum:   " << tot_str << std::endl;
            throw;
        }
        std::cout << ".";

        //
        // ===== W R I T E   T A B L E =====
        //

        // Build binary table
        GFitsBinTable table = GFitsBinTable(nrows);
        table.append_column(col);
        table.insert_column(0, col10);

        // Append to file
        fits.append(&table);

        // Save FITS file
        fits.save();
        std::cout << ".";

    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to build GFitsTableBoolCol tables."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";


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
        GFitsTableBoolCol* col =
            (GFitsTableBoolCol*)((GFitsTable*)fits.hdu(1))->column("LOGICAL");

        // Check table (operator access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i)
            tot_int += (*col)(i);
        if (tot_int != sum) {
            std::cout << std::endl << "TEST ERROR: GFitsTableBoolCol::operator()";
            std::cout << std::endl << "  Reference sum: " << sum;
            std::cout << std::endl << "  Derived sum:   " << tot_int << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i)
            tot_dbl += col->real(i);
        if (!dequal(tot_dbl, double(sum))) {
            std::cout << std::endl << "TEST ERROR: GFitsTableBoolCol::real()";
            std::cout << std::endl << "  Reference sum: " << sum;
            std::cout << std::endl << "  Derived sum:   " << tot_dbl << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i)
            tot_int += col->integer(i);
        if (tot_int != sum_int) {
            std::cout << std::endl << "TEST ERROR: GFitsTableBoolCol::integer()";
            std::cout << std::endl << "  Reference sum: " << sum_int;
            std::cout << std::endl << "  Derived sum:   " << tot_int << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i)
            tot_str += ":"+col->string(i);
        if (tot_str != sum_str) {
            std::cout << std::endl << "TEST ERROR: GFitsTableBoolCol::string()";
            std::cout << std::endl << "  Reference sum: " << sum_str;
            std::cout << std::endl << "  Derived sum:   " << tot_str << std::endl;
            throw;
        }
        std::cout << ".";

        //
        // ===== L O G I C A L  1 0 =====
        //

        // Get column
        GFitsTableBoolCol* col10 =
            (GFitsTableBoolCol*)((GFitsTable*)fits.hdu(1))->column("LOGICAL10");

        // Check table (operator access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_int += (*col10)(i,j);
        }
        if (tot_int != sum10) {
            std::cout << std::endl << "TEST ERROR: GFitsTableBoolCol::operator() - 10";
            std::cout << std::endl << "  Reference sum: " << sum10;
            std::cout << std::endl << "  Derived sum:   " << tot_int << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_dbl += col10->real(i,j);
        }
        if (!dequal(tot_dbl, sum10)) {
            std::cout << std::endl << "TEST ERROR: GFitsTableBoolCol::real() - 10";
            std::cout << std::endl << "  Reference sum: " << sum10;
            std::cout << std::endl << "  Derived sum:   " << tot_dbl << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_int += col10->integer(i,j);
        }
        if (tot_int != sum10_int) {
            std::cout << std::endl << "TEST ERROR: GFitsTableBoolCol::integer() - 10";
            std::cout << std::endl << "  Reference sum: " << sum10_int;
            std::cout << std::endl << "  Derived sum:   " << tot_int << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_str += ":"+col10->string(i,j);
        }
        if (tot_str != sum10_str) {
            std::cout << std::endl << "TEST ERROR: GFitsTableBoolCol::string() - 10";
            std::cout << std::endl << "  Reference sum: " << sum10_str;
            std::cout << std::endl << "  Derived sum:   " << tot_str << std::endl;
            throw;
        }
        std::cout << ".";
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to read GFitsTableBoolCol tables."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ". ok." << std::endl;

}


/***************************************************************************
 * @brief Test bit FITS binary table
 ***************************************************************************/
void test_bintable_bit(void)
{
    // Set filename
    std::string filename = "test_bintable_bit.fits";

    // Dump header
    std::cout << "Test GFitsTableBitCol: ";

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
        // ===== B I T =====
        //

        // Set table
        GFitsTableBitCol col = GFitsTableBitCol("BIT", nrows);
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
        std::cout << ".";

        // Check table (operator access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i)
            tot_int += col(i);
        if (tot_int != sum) {
            std::cout << std::endl << "TEST ERROR: GFitsTableBitCol::operator()";
            std::cout << std::endl << "  Reference sum: " << sum;
            std::cout << std::endl << "  Derived sum:   " << tot_int << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i)
            tot_dbl += col.real(i);
        if (!dequal(tot_dbl, double(sum))) {
            std::cout << std::endl << "TEST ERROR: GFitsTableBitCol::real()";
            std::cout << std::endl << "  Reference sum: " << sum;
            std::cout << std::endl << "  Derived sum:   " << tot_dbl << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i)
            tot_int += col.integer(i);
        if (tot_int != sum_int) {
            std::cout << std::endl << "TEST ERROR: GFitsTableBitCol::integer()";
            std::cout << std::endl << "  Reference sum: " << sum_int;
            std::cout << std::endl << "  Derived sum:   " << tot_int << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i)
            tot_str += ":"+col.string(i);
        if (tot_str != sum_str) {
            std::cout << std::endl << "TEST ERROR: GFitsTableBitCol::string()";
            std::cout << std::endl << "  Reference sum: " << sum_str;
            std::cout << std::endl << "  Derived sum:   " << tot_str << std::endl;
            throw;
        }
        std::cout << ".";

        //
        // ===== B I T  1 0 =====
        //

        // Set table
        GFitsTableBitCol col10 = GFitsTableBitCol("BIT10", nrows, nvec);
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
        std::cout << ".";

        // Check table (operator access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_int += col10(i,j);
        }
        if (tot_int != sum10) {
            std::cout << std::endl << "TEST ERROR: GFitsTableBitCol::operator() - 10";
            std::cout << std::endl << "  Reference sum: " << sum10;
            std::cout << std::endl << "  Derived sum:   " << tot_int << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_dbl += col10.real(i,j);
        }
        if (!dequal(tot_dbl, sum10)) {
            std::cout << std::endl << "TEST ERROR: GFitsTableBitCol::real() - 10";
            std::cout << std::endl << "  Reference sum: " << sum10;
            std::cout << std::endl << "  Derived sum:   " << tot_dbl << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_int += col10.integer(i,j);
        }
        if (tot_int != sum10_int) {
            std::cout << std::endl << "TEST ERROR: GFitsTableBitCol::integer() - 10";
            std::cout << std::endl << "  Reference sum: " << sum10_int;
            std::cout << std::endl << "  Derived sum:   " << tot_int << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_str += ":"+col10.string(i,j);
        }
        if (tot_str != sum10_str) {
            std::cout << std::endl << "TEST ERROR: GFitsTableBitCol::string() - 10";
            std::cout << std::endl << "  Reference sum: " << sum10_str;
            std::cout << std::endl << "  Derived sum:   " << tot_str << std::endl;
            throw;
        }
        std::cout << ".";

        //
        // ===== W R I T E   T A B L E =====
        //

        // Build binary table
        GFitsBinTable table = GFitsBinTable(nrows);
        table.append_column(col);
        table.insert_column(0, col10);

        // Append to file
        fits.append(&table);

        // Save FITS file
        fits.save();
        std::cout << ".";

    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to build GFitsTableBitCol tables."
                  << std::endl;
        std::cout << e.what() << std::endl;
        throw;
    }
    std::cout << ".";


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
        // ===== B I T =====
        //

        // Get column
        GFitsTableBitCol* col =
            (GFitsTableBitCol*)((GFitsTable*)fits.hdu(1))->column("BIT");

        // Check table (operator access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i)
            tot_int += (*col)(i);
        if (tot_int != sum) {
            std::cout << std::endl << "TEST ERROR: GFitsTableBitCol::operator()";
            std::cout << std::endl << "  Reference sum: " << sum;
            std::cout << std::endl << "  Derived sum:   " << tot_int << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i)
            tot_dbl += col->real(i);
        if (!dequal(tot_dbl, double(sum))) {
            std::cout << std::endl << "TEST ERROR: GFitsTableBitCol::real()";
            std::cout << std::endl << "  Reference sum: " << sum;
            std::cout << std::endl << "  Derived sum:   " << tot_dbl << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i)
            tot_int += col->integer(i);
        if (tot_int != sum_int) {
            std::cout << std::endl << "TEST ERROR: GFitsTableBitCol::integer()";
            std::cout << std::endl << "  Reference sum: " << sum_int;
            std::cout << std::endl << "  Derived sum:   " << tot_int << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i)
            tot_str += ":"+col->string(i);
        if (tot_str != sum_str) {
            std::cout << std::endl << "TEST ERROR: GFitsTableBitCol::string()";
            std::cout << std::endl << "  Reference sum: " << sum_str;
            std::cout << std::endl << "  Derived sum:   " << tot_str << std::endl;
            throw;
        }
        std::cout << ".";

        //
        // ===== B I T  1 0 =====
        //

        // Get column
        GFitsTableBitCol* col10 =
            (GFitsTableBitCol*)((GFitsTable*)fits.hdu(1))->column("BIT10");

        // Check table (operator access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_int += (*col10)(i,j);
        }
        if (tot_int != sum10) {
            std::cout << std::endl << "TEST ERROR: GFitsTableBitCol::operator() - 10";
            std::cout << std::endl << "  Reference sum: " << sum10;
            std::cout << std::endl << "  Derived sum:   " << tot_int << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (real access)
        tot_dbl = 0.0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_dbl += col10->real(i,j);
        }
        if (!dequal(tot_dbl, sum10)) {
            std::cout << std::endl << "TEST ERROR: GFitsTableBitCol::real() - 10";
            std::cout << std::endl << "  Reference sum: " << sum10;
            std::cout << std::endl << "  Derived sum:   " << tot_dbl << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (int access)
        tot_int = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_int += col10->integer(i,j);
        }
        if (tot_int != sum10_int) {
            std::cout << std::endl << "TEST ERROR: GFitsTableBitCol::integer() - 10";
            std::cout << std::endl << "  Reference sum: " << sum10_int;
            std::cout << std::endl << "  Derived sum:   " << tot_int << std::endl;
            throw;
        }
        std::cout << ".";

        // Check table (string access)
        tot_str.clear();
        for (int i = 0; i < nrows; ++i) {
            for (int j = 0; j < nvec; ++j)
                tot_str += ":"+col10->string(i,j);
        }
        if (tot_str != sum10_str) {
            std::cout << std::endl << "TEST ERROR: GFitsTableBitCol::string() - 10";
            std::cout << std::endl << "  Reference sum: " << sum10_str;
            std::cout << std::endl << "  Derived sum:   " << tot_str << std::endl;
            throw;
        }
        std::cout << ".";
    }
    catch (std::exception &e) {
        std::cout << std::endl
                  << "TEST ERROR: Unable to read GFitsTableBitCol tables."
                  << std::endl;
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
    test_create();
    test_image_double();
    test_bintable_bit();
    test_bintable_logical();
    test_bintable_string();
    test_bintable_double();
    test_bintable_float();
    test_bintable_ushort();
    test_bintable_short();
    test_bintable_ulong();
    test_bintable_long();
    test_bintable_longlong();
    
    // Return
    return 0;
}
