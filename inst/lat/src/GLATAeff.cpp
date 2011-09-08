/***************************************************************************
 *                 GLATAeff.cpp  -  Fermi LAT effective area               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2011 by Jurgen Knodlseder                           *
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
 * @file GLATAeff.cpp
 * @brief Fermi LAT effective area class implementation
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GLATAeff.hpp"
#include "GLATException.hpp"
#include "GTools.hpp"
#include "GFitsBinTable.hpp"
#include "GFitsTableCol.hpp"
#include "GFitsTableFloatCol.hpp"
#include "GException.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_READ                            "GLATAeff::read(const GFits* file)"
#define G_READ_AEFF                        "GLATAeff::read_aeff(GFitsTable*)"

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
GLATAeff::GLATAeff(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief File constructor
 *
 * @param[in] filename FITS file name.
 *
 * Construct instance by loading the effective area information from FITS
 * file.
 ***************************************************************************/
GLATAeff::GLATAeff(const std::string filename)
{
    // Initialise class members
    init_members();

    // Load effective area from file
    load(filename);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] aeff Effective area.
 ***************************************************************************/
GLATAeff::GLATAeff(const GLATAeff& aeff)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(aeff);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GLATAeff::~GLATAeff(void)
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
 * @param[in] aeff Effective area.
 ***************************************************************************/
GLATAeff& GLATAeff::operator= (const GLATAeff& aeff)
{
    // Execute only if object is not identical
    if (this != &aeff) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(aeff);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/***********************************************************************//**
 * @brief Return effective area in units of cm2
 *
 * @param[in] logE Log10 of the true photon energy (MeV).
 * @param[in] ctheta Cosine of zenith angle.
 ***************************************************************************/
double GLATAeff::operator() (const double& logE, const double& ctheta)
{
    // Get effective area value
    double aeff = (ctheta >= m_min_ctheta) 
                  ? m_aeff_bins.interpolate(logE, ctheta, m_aeff) : 0.0;

    // Make sure that effective area is not negative
    if (aeff < 0.0) {
        aeff = 0.0;
    }
    
    // Return effective area value
    return aeff;
}


/***********************************************************************//**
 * @brief Return effective area in units of cm2
 *
 * @param[in] logE Log10 of the true photon energy (MeV).
 * @param[in] ctheta Cosine of zenith angle.
 * @param[in] phi Azimuth angle.
 *
 * @todo Phi integration not yet implemented.
 ***************************************************************************/
double GLATAeff::operator() (const double& logE, const double& ctheta,
                             const double& phi)
{
    // Get effective area value
    double aeff = (ctheta >= m_min_ctheta) 
                  ? m_aeff_bins.interpolate(logE, ctheta, m_aeff) : 0.0;

    // Make sure that effective area is not negative
    if (aeff < 0.0) {
        aeff = 0.0;
    }

    // Return effective area value
    return aeff;
}


/***********************************************************************//**
 * @brief Return effective area (units: cm2).
 *
 * @param[in] srcDir True photon direction.
 * @param[in] srcEng True energy of photon.
 * @param[in] srcTime True photon arrival time.
 * @param[in] pnt Instrument pointing.
 *
 * @todo Not yet implemented.
 ***************************************************************************/
double GLATAeff::operator() (const GSkyDir& srcDir, const GEnergy& srcEng,
                             const GTime& srcTime, const GLATPointing& pnt)
{
    // Return Aeff value
    return 1.0;
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
void GLATAeff::clear(void)
{
    // Free class members
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
***************************************************************************/
GLATAeff* GLATAeff::clone(void) const
{
    return new GLATAeff(*this);
}


/***********************************************************************//**
 * @brief Load effective area from FITS file
 *
 * @param[in] filename FITS file.
 ***************************************************************************/
void GLATAeff::load(const std::string filename)
{
    // Open FITS file
    GFits fits(filename);

    // Read effective area from file
    read(&fits);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save effective area into FITS file
 *
 * @param[in] filename FITS file.
 * @param[in] clobber Overwrite existing file?.
 ***************************************************************************/
void GLATAeff::save(const std::string filename, bool clobber)
{
    // Open FITS file
    GFits fits(filename, true);

    // Write effective area into file
    write(fits);

    // Close FITS file
    fits.save(clobber);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read effective area from FITS file
 *
 * @param[in] fits FITS file pointer.
 *
 * @exception GException::fits_hdu_not_found
 *            Effective area HDU not found in FITS file
 ***************************************************************************/
void GLATAeff::read(const GFits* fits)
{
    // Clear instance
    clear();

    // Get pointer to effective area HDU
    GFitsTable* hdu_aeff = fits->table("EFFECTIVE AREA");
    if (hdu_aeff == NULL)
        throw GException::fits_hdu_not_found(G_READ, "EFFECTIVE AREA");

    // Read effective area
    read_aeff(hdu_aeff);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write effective area into FITS file
 *
 * @param[in] fits FITS file.
 *
 * @todo Write also phi and rate HDUs if they exist.
 ***************************************************************************/
void GLATAeff::write(GFits& fits) const
{
    // Write effective area
    write_aeff(fits);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set minimum cos(theta) angle for effective area access
 *
 * @param[in] ctheta Cosine of maximum zenith angle.
 ***************************************************************************/
void GLATAeff::costhetamin(const double& ctheta)
{
    // Set minimum cos(theta) value
    m_min_ctheta = ctheta;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print effective area information
 ***************************************************************************/
std::string GLATAeff::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GLATAeff ===");
    result.append("\n"+parformat("Number of energy bins")+str(nenergies()));
    result.append("\n"+parformat("Number of cos theta bins")+str(ncostheta()));

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                            Private methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GLATAeff::init_members(void)
{
    // Initialise members
    m_aeff_bins.clear();
    m_aeff.clear();
    m_min_ctheta = 0.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] aeff Effective area.
 ***************************************************************************/
void GLATAeff::copy_members(const GLATAeff& aeff)
{
    // Copy attributes
    m_aeff_bins  = aeff.m_aeff_bins;
    m_aeff       = aeff.m_aeff;
    m_min_ctheta = aeff.m_min_ctheta;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GLATAeff::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Read effective area from FITS table
 *
 * @param[in] hdu FITS table pointer.
 *
 * @exception GLATException::inconsistent_response
 *            Inconsistent response table encountered
 *
 * The effective area is converted into units of cm2.
 ***************************************************************************/
void GLATAeff::read_aeff(const GFitsTable* hdu)
{
    // Clear array
    m_aeff.clear();

    // Get energy and cos theta bins in response table
    m_aeff_bins.read(hdu);

    // Set minimum cos(theta)
    m_min_ctheta = m_aeff_bins.costheta_lo(0);

    // Continue only if there are effective area bins
    int size = m_aeff_bins.size();
    if (size > 0) {

        // Allocate arrays
        m_aeff.reserve(size);

        // Get pointer to effective area column
        const GFitsTableCol* ptr = &(*hdu)["EFFAREA"];

        // Check consistency of effective area table
        int num = ptr->number();
        if (num != size) {
            throw GLATException::inconsistent_response(G_READ_AEFF, num, size);
        }

        // Copy data and convert from m2 into cm2
        for (int i = 0; i < size; ++i) {
            m_aeff.push_back(ptr->real(0,i) * 1.0e4);
        }

    } // endif: there were effective area bins

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write effective area into FITS file
 *
 * @param[in] file FITS file.
 *
 * This method does not write anything if the instance is empty.
 ***************************************************************************/
void GLATAeff::write_aeff(GFits& file) const
{
    // Continue only if there are bins
    int size = m_aeff_bins.size();
    if (size > 0) {

        // Create new binary table
        GFitsBinTable* hdu_aeff = new GFitsBinTable;

        // Set table attributes
        hdu_aeff->extname("EFFECTIVE AREA");

        // Write boundaries into table
        m_aeff_bins.write(hdu_aeff);

        // Allocate floating point vector columns
        GFitsTableFloatCol col_aeff = GFitsTableFloatCol("EFFAREA",  1, size);

        // Fill columns
        for (int i = 0; i < size; ++i)
            col_aeff(0,i) = m_aeff[i];

        // Append columns to table
        hdu_aeff->append_column(col_aeff);

        // Append HDU to FITS file
        file.append(*hdu_aeff);

        // Free binary table
        delete hdu_aeff;

    } // endif: there were data to write

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Output operator
 *
 * @param[in] os Output stream.
 * @param[in] aeff Effective area.
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GLATAeff& aeff)
{
     // Write effective area in output stream
    os << aeff.print();

    // Return output stream
    return os;
}


/***********************************************************************//**
 * @brief Log operator
 *
 * @param[in] log Logger.
 * @param[in] aeff Effective area.
 ***************************************************************************/
GLog& operator<< (GLog& log, const GLATAeff& aeff)
{
    // Write effective area into logger
    log << aeff.print();

    // Return logger
    return log;
}
