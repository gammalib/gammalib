/***************************************************************************
 *       GLATPsfV3.cpp  -  Fermi/LAT point spread function version 3       *
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
 * @file GLATPsfV3.cpp
 * @brief Fermi/LAT point spread function version 3 class implementation
 * @author J. Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GLATPsfV3.hpp"
#include "GLATException.hpp"
#include "GTools.hpp"
#include "GFitsBinTable.hpp"
#include "GFitsTableFloatCol.hpp"
#include "GIntegral.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_READ                                 "GLATPsfV3::read(GFitsTable*)"

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
GLATPsfV3::GLATPsfV3(void)
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
GLATPsfV3::GLATPsfV3(const GLATPsfV3& psf)
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
GLATPsfV3::~GLATPsfV3(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Operators                                =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] psf Point spread function.
 ***************************************************************************/
GLATPsfV3& GLATPsfV3::operator= (const GLATPsfV3& psf)
{
    // Execute only if object is not identical
    if (this != &psf) {

        // Copy base class members
        this->GLATPsfBase::operator=(psf);

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

/***********************************************************************//**
 * @brief Clear instance
 *
 * This method properly resets the object to an initial state.
 ***************************************************************************/
void GLATPsfV3::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GLATPsfBase::free_members();

    // Initialise members
    this->GLATPsfBase::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
***************************************************************************/
GLATPsfV3* GLATPsfV3::clone(void) const
{
    return new GLATPsfV3(*this);
}


/***********************************************************************//**
 * @brief Read point spread function from FITS table
 *
 * @param[in] hdu FITS table pointer.
 *
 * @exception GLATException::inconsistent_response
 *            Inconsistent response table encountered
 *
 * Reads point spread function information from FITS HDU. In addition to the
 * energy and costheta binning information, 6 columns are expected:
 * NCORE, NTAIL, SCORE, STAIL, GCORE, and GTAIL.
 ***************************************************************************/
void GLATPsfV3::read(const GFitsTable* hdu)
{
    // Clear arrays
    m_ncore.clear();
    m_ntail.clear();
    m_score.clear();
    m_stail.clear();
    m_gcore.clear();
    m_gtail.clear();

    // Get energy and cos theta binning
    m_rpsf_bins.read(hdu);

    // Set minimum cos(theta)
    m_min_ctheta = m_rpsf_bins.costheta_lo(0);

    // Continue only if there are bins
    int size = m_rpsf_bins.size();
    if (size > 0) {

        // Allocate arrays
        m_ncore.reserve(size);
        m_ntail.reserve(size);
        m_score.reserve(size);
        m_stail.reserve(size);
        m_gcore.reserve(size);
        m_gtail.reserve(size);

        // Get pointer to columns
        const GFitsTableCol* ncore = &(*hdu)["NCORE"];
        const GFitsTableCol* ntail = &(*hdu)["NTAIL"];
        const GFitsTableCol* score = &(*hdu)["SCORE"];
        const GFitsTableCol* stail = &(*hdu)["STAIL"];
        const GFitsTableCol* gcore = &(*hdu)["GCORE"];
        const GFitsTableCol* gtail = &(*hdu)["GTAIL"];

        // Check consistency of columns
        if (ncore->number() != size) {
            throw GLATException::inconsistent_response(G_READ,
                                                       ncore->number(), size);
        }
        if (ntail->number() != size) {
            throw GLATException::inconsistent_response(G_READ,
                                                       ntail->number(), size);
        }
        if (score->number() != size) {
            throw GLATException::inconsistent_response(G_READ,
                                                       score->number(), size);
        }
        if (stail->number() != size) {
            throw GLATException::inconsistent_response(G_READ,
                                                       stail->number(), size);
        }
        if (gcore->number() != size) {
            throw GLATException::inconsistent_response(G_READ,
                                                       gcore->number(), size);
        }
        if (gtail->number() != size) {
            throw GLATException::inconsistent_response(G_READ,
                                                       gtail->number(), size);
        }

        // Copy data
        for (int i = 0; i < size; ++i) {
            m_ncore.push_back(ncore->real(0,i));
            m_ntail.push_back(ntail->real(0,i));
            m_score.push_back(score->real(0,i));
            m_stail.push_back(stail->real(0,i));
            m_gcore.push_back(gcore->real(0,i));
            m_gtail.push_back(gtail->real(0,i));
        }

    } // endif: there were bins

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write point spread function into FITS file
 *
 * @param[in] file FITS file.
 *
 * Writes the PSF into the extension "RPSF" of a FITS file. This method
 * does not check if a "RPSF" extension exists so far, it simply adds one
 * each time it is called.
 *
 * Nothing is done if the PSF size is 0.
 *
 * @todo Check if a RPSF extension exists already in FITS file
 ***************************************************************************/
void GLATPsfV3::write(GFits& file) const
{
    // Continue only if there are bins
    int size = m_rpsf_bins.size();
    if (size > 0) {

        // Create new binary table
        GFitsBinTable* hdu_rpsf = new GFitsBinTable;

        // Set table attributes
        hdu_rpsf->extname("RPSF");

        // Write boundaries into table
        m_rpsf_bins.write(hdu_rpsf);

        // Allocate floating point vector columns
        GFitsTableFloatCol col_ncore = GFitsTableFloatCol("NCORE",  1, size);
        GFitsTableFloatCol col_ntail = GFitsTableFloatCol("NTAIL",  1, size);
        GFitsTableFloatCol col_score = GFitsTableFloatCol("SCORE",  1, size);
        GFitsTableFloatCol col_stail = GFitsTableFloatCol("STAIL",  1, size);
        GFitsTableFloatCol col_gcore = GFitsTableFloatCol("GCORE",  1, size);
        GFitsTableFloatCol col_gtail = GFitsTableFloatCol("GTAIL",  1, size);

        // Fill columns
        for (int i = 0; i < size; ++i) {
            col_ncore(0,i) = m_ncore[i];
            col_ntail(0,i) = m_ntail[i];
            col_score(0,i) = m_score[i];
            col_stail(0,i) = m_stail[i];
            col_gcore(0,i) = m_gcore[i];
            col_gtail(0,i) = m_gtail[i];
        }

        // Append columns to table
        hdu_rpsf->append_column(col_ncore);
        hdu_rpsf->append_column(col_ntail);
        hdu_rpsf->append_column(col_score);
        hdu_rpsf->append_column(col_stail);
        hdu_rpsf->append_column(col_gcore);
        hdu_rpsf->append_column(col_gtail);

        // Append HDU to FITS file
        file.append(*hdu_rpsf);

        // Free binary table
        delete hdu_rpsf;

    } // endif: there were data to write

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return point spread function value
 *
 * @param[in] offset Offset angle (deg).
 * @param[in] logE Log10 of the true photon energy (MeV).
 * @param[in] ctheta Cosine of zenith angle.
 *
 * @todo Implement method.
 ***************************************************************************/
double GLATPsfV3::psf(const double& offset, const double& logE,
                            const double& ctheta)
{
    // Initialise response
    double psf = 0.0;

    // Return point spread function
    return psf;
}


/*==========================================================================
 =                                                                         =
 =                            Private methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GLATPsfV3::init_members(void)
{
    // Initialise members
    m_ncore.clear();
    m_ntail.clear();
    m_score.clear();
    m_stail.clear();
    m_gcore.clear();
    m_gtail.clear();
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] psf Point spread function.
 ***************************************************************************/
void GLATPsfV3::copy_members(const GLATPsfV3& psf)
{
    // Copy members
    m_ncore = psf.m_ncore;
    m_ntail = psf.m_ntail;
    m_score = psf.m_score;
    m_stail = psf.m_stail;
    m_gcore = psf.m_gcore;
    m_gtail = psf.m_gtail;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GLATPsfV3::free_members(void)
{
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/
