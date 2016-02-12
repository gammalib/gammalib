/***************************************************************************
 *                GCTAAeff2D.cpp - CTA 2D effective area class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2016 by Juergen Knoedlseder                         *
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
 * @file GCTAAeff2D.hpp
 * @brief CTA 2D effective area class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GException.hpp"
#include "GFilename.hpp"
#include "GFits.hpp"
#include "GFitsBinTable.hpp"
#include "GFitsTable.hpp"
#include "GCTAAeff2D.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_READ                                "GCTAAeff2D::read(GFitsTable&)"

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
 *
 * Constructs empty effective area.
 ***************************************************************************/
GCTAAeff2D::GCTAAeff2D(void) : GCTAAeff()
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
 * Constructs effective area from a FITS file.
 ***************************************************************************/
GCTAAeff2D::GCTAAeff2D(const std::string& filename) : GCTAAeff()
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
 *
 * Constructs effective area by copying from another effective area.
 ***************************************************************************/
GCTAAeff2D::GCTAAeff2D(const GCTAAeff2D& aeff) : GCTAAeff(aeff)
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
 *
 * Destructs effective area.
 ***************************************************************************/
GCTAAeff2D::~GCTAAeff2D(void)
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
 * @return Effective area.
 *
 * Assigns effective area.
 ***************************************************************************/
GCTAAeff2D& GCTAAeff2D::operator=(const GCTAAeff2D& aeff)
{
    // Execute only if object is not identical
    if (this != &aeff) {

        // Copy base class members
        this->GCTAAeff::operator=(aeff);

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
 * @param[in] logE Log10 of the true photon energy (TeV).
 * @param[in] theta Offset angle in camera system (rad) (default: 0.0).
 * @param[in] phi Azimuth angle in camera system (rad).
 * @param[in] zenith Zenith angle in Earth system (rad).
 * @param[in] azimuth Azimuth angle in Earth system (rad).
 * @param[in] etrue Use true energy (default: true).
 * @return Effective area in cm2.
 *
 * Returns the effective area in units of cm2 for a given energy and
 * offset angle. The effective area is bi-linearly interpolated in the
 * log10(energy) - offset angle plane. The method assures that the effective
 * area value never becomes negative.
 *
 * The method supports true and reconstructed energies for logE. To access
 * the effective area as function of true energy, specify etrue=true
 * (this is the default). The obtained the effective area as function of
 * reconstructed energy, specify etrue=false.
 ***************************************************************************/
double GCTAAeff2D::operator()(const double& logE, 
                              const double& theta, 
                              const double& phi,
                              const double& zenith,
                              const double& azimuth,
                              const bool&   etrue) const
{
    // Set parameter index
    int index = (etrue) ? m_inx_aeff : m_inx_aeff_reco;

    // Get effective area value in cm2
    double aeff = (m_inx_energy == 0) ? m_aeff(index, logE, theta) :
                                        m_aeff(index, theta, logE);

    // Make sure that effective area is not negative
    if (aeff < 0.0) {
        aeff = 0.0;
    }
    
    // Return effective area value
    return aeff;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear effective area.
 *
 * Clears effective area.
 ***************************************************************************/
void GCTAAeff2D::clear(void)
{
    // Free class members
    free_members();
    this->GCTAAeff::free_members();

    // Initialise members
    this->GCTAAeff::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone effective area
 *
 * @return Deep copy of effective area.
 *
 * Returns a pointer to a deep copy of the effective area.
 ***************************************************************************/
GCTAAeff2D* GCTAAeff2D::clone(void) const
{
    return new GCTAAeff2D(*this);
}


/***********************************************************************//**
 * @brief Read effective area from FITS table
 *
 * @param[in] table FITS table.
 *
 * @exception GException::invalid_value
 *            Response table is not two-dimensional.
 *
 * Reads the effective area form the FITS @p table. The following column
 * names are mandatory:
 *
 *     ENERG_LO - Energy lower bin boundaries
 *     ENERG_HI - Energy upper bin boundaries
 *     THETA_LO - Offset angle lower bin boundaries
 *     THETA_HI - Offset angle upper bin boundaries
 *     EFFAREA  - Effective area
 *
 * In addition, the following column names are optional:
 *
 *     EFFAREA_RECO - Effective area as function of reconstructed energy
 *
 * The data are stored in the m_aeff member. The energy axis will be set to
 * log10, the offset angle axis to radians.
 *
 * @todo Analyse the unit of the parameter axis to determine the conversion
 * factor for the effective areas. For the moment they are hard wired. 
 ***************************************************************************/
void GCTAAeff2D::read(const GFitsTable& table)
{
    // Clear response table
    m_aeff.clear();

    // Read effective area table
    m_aeff.read(table);

    // Get mandatory indices (throw exception if not found)
    m_inx_energy = m_aeff.axis("ENERG");
    m_inx_theta  = m_aeff.axis("THETA");
    m_inx_aeff   = m_aeff.table("EFFAREA");

    // Get optional index (use "EFFAREA" if "EFFAREA_RECO" does not exist)
    m_inx_aeff_reco = (m_aeff.has_table("EFFAREA_RECO")) ?
                       m_aeff.table("EFFAREA_RECO") : m_inx_aeff;

    // Throw an exception if the table is not two-dimensional
    if (m_aeff.axes() != 2) {
        std::string msg = "Expected two-dimensional effective area response "
                          "table but found "+gammalib::str(m_aeff.axes())+
                          " dimensions. Please specify a two-dimensional "
                          "effective area.";
        throw GException::invalid_value(G_READ, msg);
    }

    // Set energy axis to logarithmic scale
    m_aeff.axis_log10(m_inx_energy);

    // Set offset angle axis to radians
    m_aeff.axis_radians(m_inx_theta);

    // Convert effective areas from m2 to cm2
    for (int i = 0; i < m_aeff.tables(); ++i) {
        m_aeff.scale(i, 1.0e4);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write effective area into FITS binary table
 *
 * @param[in] table FITS binary table.
 *
 * Writes effective area into a FITS binary @p table.
 *
 * @todo Add keywords.
 ***************************************************************************/
void GCTAAeff2D::write(GFitsBinTable& table) const
{
    // Create a copy of the response table
    GCTAResponseTable aeff(m_aeff);

    // Convert area from cm2 to m2
    for (int i = 0; i < aeff.tables(); ++i) {
        aeff.scale(i, 1.0e-4);
    }

    // Write response table
    aeff.write(table);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Load effective area from FITS file
 *
 * @param[in] filename FITS file name.
 *
 * Loads the effective area from a FITS file.
 *
 * If no extension name is provided, the effective area will be loaded from
 * the "EFFECTIVE AREA" extension.
 ***************************************************************************/
void GCTAAeff2D::load(const std::string& filename)
{
    // Create file name
    GFilename fname(filename);

    // Allocate FITS file
    GFits file;

    // Open FITS file
    file.open(fname.filename());

    // Get effective area table
    const GFitsTable& table = *file.table(fname.extname("EFFECTIVE AREA"));

    // Read effective area from table
    read(table);

    // Close FITS file
    file.close();

    // Store filename
    m_filename = filename;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save effectiva area into FITS file
 *
 * @param[in] filename FITS file name.
 * @param[in] clobber Overwrite existing file? (default: false)
 *
 * Save the effectiva area into a FITS file.
 *
 * If no extension name is provided, the effective area will be saved into
 * the "EFFECTIVE AREA" extension.
 ***************************************************************************/
void GCTAAeff2D::save(const std::string& filename, const bool& clobber) const
{
    // Create file name
    GFilename fname(filename);

    // Create binary table
    GFitsBinTable table;
    table.extname(fname.extname("EFFECTIVE AREA"));

    // Write the Effective area table
    write(table);

    // Create FITS file, append table, and write into the file
    GFits fits;
    fits.append(table);
    fits.saveto(fname.filename(), clobber);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return maximum effective area at a given energy in cm2
 *
 * @param[in] logE Log10 of the true photon energy (TeV).
 * @param[in] zenith Zenith angle in Earth system (rad).
 * @param[in] azimuth Azimuth angle in Earth system (rad).
 * @param[in] etrue Use true energy (default: true).
 * @return Maximum effective area (cm2).
 *
 * Returns the maximum effective area for a given energy, zenith and azimuth
 * angle in units of cm2.
 ***************************************************************************/
double GCTAAeff2D::max(const double& logE,
                       const double& zenith,
                       const double& azimuth,
                       const bool&   etrue) const
{
    // Initialise maximum effective area
    double max_aeff = 0.0;

    // Set parameter index
    int index = (etrue) ? m_inx_aeff : m_inx_aeff_reco;

    // Get number of theta bins
    int n_theta = m_aeff.axis_bins(m_inx_theta);

    // Loop over theta values
    for (int i = 0; i < n_theta; ++i) {

        // Compute lower and upper theta bin values
        double theta_lo = m_aeff.axis_lo(m_inx_theta, i);
        double theta_hi = m_aeff.axis_hi(m_inx_theta, i);

        // Get effective area value in cm2
        double aeff_lo = (m_inx_energy == 0) ? m_aeff(index, logE, theta_lo) :
                                               m_aeff(index, theta_lo, logE);
        double aeff_hi = (m_inx_energy == 0) ? m_aeff(index, logE, theta_hi) :
                                               m_aeff(index, theta_hi, logE);

        // Update maximum effective area if larger than current maximum
        // effective area
        if (aeff_lo > max_aeff) {
            max_aeff = aeff_lo;
        }
        if (aeff_hi > max_aeff) {
            max_aeff = aeff_hi;
        }

    } // endfor: loop over theta values

    // Return effective area value
    return max_aeff;
}


/***********************************************************************//**
 * @brief Print effective area information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing effective area information.
 ***************************************************************************/
std::string GCTAAeff2D::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Compute energy boundaries in TeV
        double emin = m_aeff.axis_lo(m_inx_energy,0);
        double emax = m_aeff.axis_hi(m_inx_energy,
                                     m_aeff.axis_bins(m_inx_energy)-1);

        // Compute offset angle boundaries in deg
        double omin = m_aeff.axis_lo(m_inx_theta,0);
        double omax = m_aeff.axis_hi(m_inx_theta,
                                     m_aeff.axis_bins(m_inx_theta)-1);

        // Append header
        result.append("=== GCTAAeff2D ===");

        // Append information
        result.append("\n"+gammalib::parformat("Filename")+m_filename);
        result.append("\n"+gammalib::parformat("Number of energy bins") +
                      gammalib::str(m_aeff.axis_bins(m_inx_energy)));
        result.append("\n"+gammalib::parformat("Number of offset bins") +
                      gammalib::str(m_aeff.axis_bins(m_inx_theta)));
        result.append("\n"+gammalib::parformat("Log10(Energy) range"));
        result.append(gammalib::str(emin)+" - "+gammalib::str(emax)+" TeV");
        result.append("\n"+gammalib::parformat("Offset angle range"));
        result.append(gammalib::str(omin)+" - "+gammalib::str(omax)+" deg");

    } // endif: chatter was not silent

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
void GCTAAeff2D::init_members(void)
{
    // Initialise members
    m_filename.clear();
    m_aeff.clear();
    m_inx_energy    = 0;
    m_inx_theta     = 1;
    m_inx_aeff      = 0;
    m_inx_aeff_reco = 1;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] aeff Effective area.
 ***************************************************************************/
void GCTAAeff2D::copy_members(const GCTAAeff2D& aeff)
{
    // Copy members
    m_filename      = aeff.m_filename;
    m_aeff          = aeff.m_aeff;
    m_inx_energy    = aeff.m_inx_energy;
    m_inx_theta     = aeff.m_inx_theta;
    m_inx_aeff      = aeff.m_inx_aeff;
    m_inx_aeff_reco = aeff.m_inx_aeff_reco;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCTAAeff2D::free_members(void)
{
    // Return
    return;
}
