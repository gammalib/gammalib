/***************************************************************************
 *                 GCTAAeff2D.cpp - CTA 2D effective area class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2015 by Juergen Knoedlseder                         *
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
 ***************************************************************************/
GCTAAeff2D& GCTAAeff2D::operator= (const GCTAAeff2D& aeff)
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
 * @param[in] theta Offset angle in camera system (rad). Defaults to 0.0.
 * @param[in] phi Azimuth angle in camera system (rad). Not used in this method.
 * @param[in] zenith Zenith angle in Earth system (rad). Not used in this method.
 * @param[in] azimuth Azimuth angle in Earth system (rad). Not used in this method.
 * @param[in] etrue Use true energy (true/false). Defaults to true.
 *
 * Returns the effective area in units of cm2 for a given energy and
 * offset angle. The effective area is bi-linearily interpolated in the
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
    int index = (etrue) ? 0 : 1;

    // Get effective area value in cm2
    double aeff = m_aeff(index, logE, theta);

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
 * @brief Clear instance
 *
 * This method properly resets the object to an initial state.
 ***************************************************************************/
void GCTAAeff2D::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GCTAAeff::free_members();

    // Initialise members
    this->GCTAAeff::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
 *
 * @return Deep copy of effective area instance.
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
 *            FITS file format differs from expectation.
 *
 * Reads the effective area form the FITS @p table.
 *
 * The data are stored in m_aeff which is of type GCTAResponseTable. The
 * energy axis will be set to log10, the offset angle axis to radians.
 ***************************************************************************/
void GCTAAeff2D::read(const GFitsTable& table)
{
    // Clear response table
    m_aeff.clear();

    // Read effective area table
    m_aeff.read(table);

    // Check that axis names comply to format
    if (m_aeff.axis_lo_name(0) != "ENERG_LO" ||
        m_aeff.axis_hi_name(0) != "ENERG_HI") {
        std::string msg = "Effective area response table does not contain"
                          " \"ENERG_LO\" and \"ENERG_HI\" columns as the"
                          " first axis.";
        throw GException::invalid_value(G_READ, msg);
    }
    if (m_aeff.axis_lo_name(1) != "THETA_LO" ||
        m_aeff.axis_hi_name(1) != "THETA_HI") {
        std::string msg = "Effective area response table does not contain"
                          " \"THETA_LO\" and \"THETA_HI\" columns as the"
                          " second axis.";
        throw GException::invalid_value(G_READ, msg);
    }

    // Set energy axis to logarithmic scale
    m_aeff.axis_log10(0);

    // Set offset angle axis to radians
    m_aeff.axis_radians(1);

    // Convert effective areas from m2 to cm2
    m_aeff.scale(0, 1.0e4);
    m_aeff.scale(1, 1.0e4);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write effective area into FITS binary table
 *
 * @param[in] table FITS binary table.
 *
 * Writes effective area into the FITS binary @p table.
 *
 * @todo Add keywords.
 ***************************************************************************/
void GCTAAeff2D::write(GFitsBinTable& table) const
{
    // Create a copy of the response table
    GCTAResponseTable aeff(m_aeff);

    // Convert area from cm2 to m2
    aeff.scale(0, 1.0e-4);
    aeff.scale(1, 1.0e-4);

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
 * @brief Return maximum effective area at a given energy
 *
 * @param[in] logE Log10 of the true photon energy (TeV).
 * @param[in] zenith Zenith angle in Earth system (rad). Not used in this method.
 * @param[in] azimuth Azimuth angle in Earth system (rad). Not used in this method.
 * @param[in] etrue Use true energy (true/false). Not used.
 * @return Maximum effective area (cm2).
 ***************************************************************************/
double GCTAAeff2D::max(const double& logE,
                       const double& zenith,
                       const double& azimuth,
                       const bool&   etrue) const
{
    // Set parameter index
    int index = (etrue) ? 0 : 1;

    // Get number of theta bins
    int n_theta = m_aeff.axis(1);

    // Initialise maximum effective area
    double max_aeff = 0.0;

    // Compute stepsize
    double theta_step = m_aeff.axis_hi(1, n_theta - 1) / (double(n_theta) - 1.0);

    // Loop over theta values
    for (int i = 0; i < n_theta; ++i) {

        // Compute theta
        double theta = i * theta_step;

        // Get effective area value in cm2
        double aeff = m_aeff(index, logE, theta);

        // Update maximum effective area if larger than current maximum
        // effective area
        if (aeff > max_aeff) {
            max_aeff = aeff;
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
        double emin = m_aeff.axis_lo(0,0);
        double emax = m_aeff.axis_hi(0,m_aeff.axis(0)-1);

        // Compute offset angle boundaries in deg
        double omin = m_aeff.axis_lo(1,0);
        double omax = m_aeff.axis_hi(1,m_aeff.axis(1)-1);

        // Append header
        result.append("=== GCTAAeff2D ===");

        // Append information
        result.append("\n"+gammalib::parformat("Filename")+m_filename);
        result.append("\n"+gammalib::parformat("Number of energy bins") +
                      gammalib::str(m_aeff.axis(0)));
        result.append("\n"+gammalib::parformat("Number of offset bins") +
                      gammalib::str(m_aeff.axis(1)));
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
    m_filename = aeff.m_filename;
    m_aeff     = aeff.m_aeff;

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
