/***************************************************************************
 *         GLATResponseTable.cpp  -  GLAST LAT Response table class        *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2008 by Jurgen Knodlseder                   *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include <math.h>
#include "GException.hpp"
#include "GLATResponseTable.hpp"
#include "GFitsHDU.hpp"
#include "GFitsTableFltCol.hpp"

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Constants __________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                GLATResponseTable constructors/destructors               =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GLATResponseTable::GLATResponseTable()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param table Response table to be copied
 ***************************************************************************/
GLATResponseTable::GLATResponseTable(const GLATResponseTable& table)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(table);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GLATResponseTable::~GLATResponseTable()
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                       GLATResponseTable operators                       =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param table Response table to be assigned
 ***************************************************************************/
GLATResponseTable& GLATResponseTable::operator= (const GLATResponseTable& table)
{
    // Execute only if object is not identical
    if (this != &table) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(table);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                      GLATResponseTable public methods                   =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Load response table from FITS HDU
 *
 * @param[in] hdu Pointer to HDU from which the response table should be loaded
 *
 * The response table definition is assumed to be stored in 4 vector columns
 * with the names
 * ENERG_LO (energy bins lower boundary)
 * ENERG_HI (energy bins upper boundary)
 * CTHETA_LO (cos theta bins lower boundary)
 * CTHETA_HI (cos theta bins upper boundary)
 ***************************************************************************/
void GLATResponseTable::load(const GFitsHDU* hdu)
{
    // Get pointers to table columns
    GFitsTableCol* energy_lo = hdu->column("ENERG_LO");
    GFitsTableCol* energy_hi = hdu->column("ENERG_HI");
    GFitsTableCol* ctheta_lo = hdu->column("CTHETA_LO");
    GFitsTableCol* ctheta_hi = hdu->column("CTHETA_HI");
    
    // Check on existence of columns
    if (energy_lo == NULL || energy_hi == NULL ||
        ctheta_lo == NULL || ctheta_hi == NULL)
        return;

    // Free any existing memory
    free_members();
        
    // Extract number of bins
    m_energy_num = energy_lo->number();
    m_ctheta_num = ctheta_lo->number();
    
    // Allocate memory
    if (m_energy_num > 0) {
        m_energy_lo = new double[m_energy_num];
        m_energy_hi = new double[m_energy_num];
    }
    if (m_ctheta_num > 0) {
        m_ctheta_lo = new double[m_ctheta_num];
        m_ctheta_hi = new double[m_ctheta_num];
    }
    
    // Transfer data
    for (int i = 0; i < m_energy_num; ++i) {
        m_energy_lo[i] = energy_lo->real(0,i);
        m_energy_hi[i] = energy_hi->real(0,i);
    }
    for (int i = 0; i < m_ctheta_num; ++i) {
        m_ctheta_lo[i] = ctheta_lo->real(0,i);
        m_ctheta_hi[i] = ctheta_hi->real(0,i);
    }
    
    // Set energy nodes
    if (m_energy_num > 0) {
        double *logE = new double[m_energy_num];
        for (int i = 0; i < m_energy_num; ++i)
            logE[i] = 0.5 * (log10(m_energy_lo[i]) + log10(m_energy_hi[i]));
        m_energy.nodes(m_energy_num, logE);
        delete [] logE;
    }

    // Set cos theta nodes
    if (m_ctheta_num > 0) {
        double *ctheta = new double[m_ctheta_num];
        for (int i = 0; i < m_ctheta_num; ++i)
            ctheta[i] = 0.5 * (m_ctheta_lo[i] + m_ctheta_hi[i]);
        m_ctheta.nodes(m_ctheta_num, ctheta);
        delete [] ctheta;
    }
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Save response table into FITS HDU
 *
 * @param[in] hdu HDU into which the response table should be saved.
 ***************************************************************************/
void GLATResponseTable::save(GFitsHDU* hdu) const
{
    // Allocate boundary table with one row
    GFitsBinTable table(1);
    
    // Allocate floating point vector columns
    GFitsTableFltCol col_energy_lo = GFitsTableFltCol("ENERG_LO",  1, m_energy_num);
    GFitsTableFltCol col_energy_hi = GFitsTableFltCol("ENERG_HI",  1, m_energy_num);
    GFitsTableFltCol col_ctheta_lo = GFitsTableFltCol("CTHETA_LO", 1, m_ctheta_num);
    GFitsTableFltCol col_ctheta_hi = GFitsTableFltCol("CTHETA_HI", 1, m_ctheta_num);

    // Set column values
    for (int i = 0; i < m_energy_num; ++i) {
        col_energy_lo(0,i) = m_energy_lo[i];
        col_energy_hi(0,i) = m_energy_hi[i];
    }
    for (int i = 0; i < m_ctheta_num; ++i) {
        col_ctheta_lo(0,i) = m_ctheta_lo[i];
        col_ctheta_hi(0,i) = m_ctheta_hi[i];
    }

    // Append columns to boundary table
    table.append_column(col_energy_lo);
    table.append_column(col_energy_hi);
    table.append_column(col_ctheta_lo);
    table.append_column(col_ctheta_hi);

    // Set HDU
    *hdu = GFitsHDU(table);
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Return table index
 *
 * @param[in] ie Energy bin (starting from 0)
 * @param[in] ic cos theta bin (starting from 0)
 ***************************************************************************/
int GLATResponseTable::index(const int& ie, const int& ic) const
{
    return ic * m_energy_num + ie;
}


/***********************************************************************//**
 * @brief Return mean energy of bin
 *
 * @param[in] ie Index of energy bin (starting from 0)
 *
 * The mean energy is actually computed from the logarithmic average of 
 * lower and upper boundaries:
 * /f$E=10^{0.5 \times (\log{E_{\rm min}} + \log{E_{\rm max}})}/f$
 ***************************************************************************/
double GLATResponseTable::energy(const int& ie) const
{
    // Determine mean energy of bin
    double mean   = 0.5 * (log10(m_energy_lo[ie]) + log10(m_energy_hi[ie]));
    double energy = pow(10.0, mean);

    // Return mean energy
    return energy;
}


/***********************************************************************//**
 * @brief Perform bi-linear interpolation of 2D array
 *
 * @param[in] logE Base 10 logarithm of the energy (MeV)
 * @param[in] ctheta cos(theta)
 * @param[in] array Array to be interpolated
 ***************************************************************************/
double GLATResponseTable::interpolate(const double& logE, 
                                      const double& ctheta, 
                                      const double* array)
{
    // Flag no change of values
    int change = 0;
    
    // Set energy interpolation
    if (logE != m_last_energy) {
        m_energy.set_value(logE);
        m_last_energy = logE;
        change        = 1;
    }

    // Set cos(theta) interpolation
    if (ctheta != m_last_ctheta) {
        m_ctheta.set_value(ctheta);
        m_last_ctheta = ctheta;
        change        = 1;
    }
    
    // If change occured then update interpolation indices and weighting
    // factors
    if (change) {
    
        // Set array indices for bi-linear interpolation
        int inx_ctheta_left  = m_ctheta.inx_left()  * m_energy_num;
        int inx_ctheta_right = m_ctheta.inx_right() * m_energy_num;
        m_inx1 = m_energy.inx_left()  + inx_ctheta_left;
        m_inx2 = m_energy.inx_left()  + inx_ctheta_right;
        m_inx3 = m_energy.inx_right() + inx_ctheta_left;
        m_inx4 = m_energy.inx_right() + inx_ctheta_right;

        // Set weighting factors for bi-linear interpolation
        m_wgt1 = m_energy.wgt_left()  * m_ctheta.wgt_left();
        m_wgt2 = m_energy.wgt_left()  * m_ctheta.wgt_right();
        m_wgt3 = m_energy.wgt_right() * m_ctheta.wgt_left();
        m_wgt4 = m_energy.wgt_right() * m_ctheta.wgt_right();
        
    } // endif: logE or ctheta changed
    
    // Perform bi-linear interpolation
    double value = m_wgt1 * array[m_inx1] +
                   m_wgt2 * array[m_inx2] +
                   m_wgt3 * array[m_inx3] +
                   m_wgt4 * array[m_inx4];
    
    // Return bi-linear interpolated value
    return value;
}


/***********************************************************************//**
 * @brief Perform bi-linear interpolation of 3D array
 *
 * @param[in] logE Base 10 logarithm of the energy (MeV)
 * @param[in] ctheta cos(theta)
 * @param[in] array Array to be interpolated
 * @param[in] offset Offset if 3D array in 1st dimension
 * @param[in] size Size of 3D array in 1st dimension
 ***************************************************************************/
double GLATResponseTable::interpolate(const double& logE, 
                                      const double& ctheta, 
                                      const double* array,
                                      const int&    offset,
                                      const int&    size)
{
    // Flag no change of values
    int change = 0;
    
    // Set energy interpolation
    if (logE != m_last_energy) {
        m_energy.set_value(logE);
        m_last_energy = logE;
        change        = 1;
    }

    // Set cos(theta) interpolation
    if (ctheta != m_last_ctheta) {
        m_ctheta.set_value(ctheta);
        m_last_ctheta = ctheta;
        change        = 1;
    }
    
    // If change occured then update interpolation indices and weighting
    // factors
    if (change) {
    
        // Set array indices for bi-linear interpolation
        int inx_ctheta_left  = m_ctheta.inx_left()  * m_energy_num;
        int inx_ctheta_right = m_ctheta.inx_right() * m_energy_num;
        m_inx1 = m_energy.inx_left()  + inx_ctheta_left;
        m_inx2 = m_energy.inx_left()  + inx_ctheta_right;
        m_inx3 = m_energy.inx_right() + inx_ctheta_left;
        m_inx4 = m_energy.inx_right() + inx_ctheta_right;

        // Set weighting factors for bi-linear interpolation
        m_wgt1 = m_energy.wgt_left()  * m_ctheta.wgt_left();
        m_wgt2 = m_energy.wgt_left()  * m_ctheta.wgt_right();
        m_wgt3 = m_energy.wgt_right() * m_ctheta.wgt_left();
        m_wgt4 = m_energy.wgt_right() * m_ctheta.wgt_right();
        
    } // endif: logE or ctheta changed
    
    // Perform bi-linear interpolation
    double value = m_wgt1 * array[m_inx1*size+offset] +
                   m_wgt2 * array[m_inx2*size+offset] +
                   m_wgt3 * array[m_inx3*size+offset] +
                   m_wgt4 * array[m_inx4*size+offset];
    
    // Return bi-linear interpolated value
    return value;
}


/*==========================================================================
 =                                                                         =
 =                      GLATResponseTable private methods                  =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GLATResponseTable::init_members(void)
{
    // Initialise members
    m_energy_num  = 0;
    m_ctheta_num  = 0;
    m_energy_lo   = NULL;
    m_energy_hi   = NULL;
    m_ctheta_lo   = NULL;
    m_ctheta_hi   = NULL;
    m_last_energy = -1.0;
    m_last_ctheta = -1.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param table Response table to be copied
 ***************************************************************************/
void GLATResponseTable::copy_members(const GLATResponseTable& table)
{
    // Copy number of bins
    m_energy_num  = table.m_energy_num;
    m_ctheta_num  = table.m_ctheta_num;
    m_energy      = table.m_energy;
    m_ctheta      = table.m_ctheta;
    m_last_energy = table.m_last_energy;
    m_last_ctheta = table.m_last_ctheta;
    
    // Copy energy bins
    if (m_energy_num > 0) {
        m_energy_lo   = new double[m_energy_num];
        m_energy_hi   = new double[m_energy_num];
        memcpy(m_energy_lo, table.m_energy_lo, m_energy_num*sizeof(double));
        memcpy(m_energy_hi, table.m_energy_hi, m_energy_num*sizeof(double));
    }

    // Copy cos theta bins
    if (m_ctheta_num > 0) {
        m_ctheta_lo = new double[m_ctheta_num];
        m_ctheta_hi = new double[m_ctheta_num];
        memcpy(m_ctheta_lo, table.m_ctheta_lo, m_ctheta_num*sizeof(double));
        memcpy(m_ctheta_hi, table.m_ctheta_hi, m_ctheta_num*sizeof(double));
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GLATResponseTable::free_members(void)
{
    // Free memory
    if (m_energy_lo != NULL) delete [] m_energy_lo;
    if (m_energy_hi != NULL) delete [] m_energy_hi;
    if (m_ctheta_lo != NULL) delete [] m_ctheta_lo;
    if (m_ctheta_hi != NULL) delete [] m_ctheta_hi;

    // Signal that memory is free
    m_energy_lo = NULL;
    m_energy_hi = NULL;
    m_ctheta_lo = NULL;
    m_ctheta_hi = NULL;

    // Return
    return;
}
