/***************************************************************************
 * GLATPsfBase.cpp  -  Fermi/LAT point spread function abstract base class *
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
 * @file GLATPsfBase.cpp
 * @brief Fermi/LAT point spread function abstract base class implementation
 * @author J. Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GLATPsfBase.hpp"
#include "GLATException.hpp"
#include "GTools.hpp"
#include "GFitsBinTable.hpp"
#include "GFitsTableFloatCol.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_READ_SCALE                   "GLATPsfBase::read_scale(GFitsTable*)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Constants __________________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GLATPsfBase::GLATPsfBase(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] psf Point spread function.
 ***************************************************************************/
GLATPsfBase::GLATPsfBase(const GLATPsfBase& psf)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(psf);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GLATPsfBase::~GLATPsfBase(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                               Operators                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] psf Point spread function.
 ***************************************************************************/
GLATPsfBase& GLATPsfBase::operator= (const GLATPsfBase& psf)
{
    // Execute only if object is not identical
    if (this != &psf) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(psf);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GLATPsfBase::init_members(void)
{
    // Initialise members
    m_front       = true;
    m_rpsf_bins.clear();
    m_scale_par1  = 0.0;
    m_scale_par2  = 0.0;
    m_scale_index = 0.0;
    m_min_ctheta  = 0.0;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] psf Point spread function.
 ***************************************************************************/
void GLATPsfBase::copy_members(const GLATPsfBase& psf)
{
    // Copy members
    m_front       = psf.m_front;
    m_rpsf_bins   = psf.m_rpsf_bins;
    m_scale_par1  = psf.m_scale_par1;
    m_scale_par2  = psf.m_scale_par2;
    m_scale_index = psf.m_scale_index;
    m_min_ctheta  = psf.m_min_ctheta;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GLATPsfBase::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Read PSF scale factors from FITS table
 *
 * @param[in] hdu FITS table pointer.
 *
 * Reads the PSF scale factors from column "PSFSCALE" of a FITS table.
 ***************************************************************************/
void GLATPsfBase::read_scale(const GFitsTable* hdu)
{
    // Get pointer to column
    const GFitsTableCol* scale = &(*hdu)["PSFSCALE"];

    // Get scaling factors
    if (front()) {
        m_scale_par1 = scale->real(0,0);
        m_scale_par2 = scale->real(0,1);
    }
    else {
        m_scale_par1 = scale->real(0,2);
        m_scale_par2 = scale->real(0,3);
    }
    m_scale_index = scale->real(0,4);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write PSF scale factors
 *
 * @param[in] file FITS file.
 *
 * Writes the PSF scale factors in the extension "PSF_SCALING_PARAMS". This
 * method appends an extension "PSF_SCALING_PARAMS" to the FITS file,
 * irrespectively of whether the extension exists already or not.
 * The scale facors are written into the column "PSFSCALE".
 *
 * @todo Check if PSF_SCALING_PARAMS exists already in FITS file
 ***************************************************************************/
void GLATPsfBase::write_scale(GFits& file) const
{
    // Create new binary table
    GFitsBinTable* hdu_scale = new GFitsBinTable;

    // Set table attributes
    hdu_scale->extname("PSF_SCALING_PARAMS");

    // Allocate floating point vector column
    GFitsTableFloatCol col_scale = GFitsTableFloatCol("PSFSCALE",  1, 5);

    // Fill columns
    if (front()) {
        col_scale(0,0) = m_scale_par1;
        col_scale(0,1) = m_scale_par2;
        col_scale(0,2) = 0.0;
        col_scale(0,3) = 0.0;
    }
    else {
        col_scale(0,0) = 0.0;
        col_scale(0,1) = 0.0;
        col_scale(0,2) = m_scale_par1;
        col_scale(0,3) = m_scale_par2;
    }
    col_scale(0,4) = m_scale_index;

    // Append column to table
    hdu_scale->append_column(col_scale);

    // Append HDU to FITS file
    file.append(*hdu_scale);

    // Free binary table
    delete hdu_scale;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return scale factor for energy (in MeV)
 *
 * @param[in] energy Photon energy (in MeV).
 ***************************************************************************/
double GLATPsfBase::scale_factor(const double& energy) const
{
    // Compute scale factor
    double f1    = m_scale_par1 * pow(0.01*energy, m_scale_index);
    double scale = sqrt(f1*f1 + m_scale_par2*m_scale_par2);

    // Return scale factor
    return scale;
}
