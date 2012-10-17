/***************************************************************************
 *                 GCTAAeff2D.hpp - CTA 2D effective area class            *
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
 * @file GCTAAeff2D.hpp
 * @brief CTA 2D effective area class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
//#include <cmath>
#include "GTools.hpp"
#include "GFitsTable.hpp"
#include "GCTAAeff2D.hpp"

/* __ Method name definitions ____________________________________________ */

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
 * Construct instance by loading the effective area information from a FITS
 * file.
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
 * @param[in] theta Offset angle (rad). Defaults to 0.0.
 * @param[in] phi Azimuth angle (rad). Not used in this method.
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
 * @brief Load effective area from FITS file
 *
 * @param[in] filename FITS file.
 *
 * This method loads the effective area information from a FITS file.
 ***************************************************************************/
void GCTAAeff2D::load(const std::string& filename)
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
 *
 * This method saves the effective area information into a FITS file.
 ***************************************************************************/
void GCTAAeff2D::save(const std::string& filename, const bool& clobber) const
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
 * Reads the effective area form the FITS file extension "EFFECTIVE AREA".
 * The data are stored in m_aeff which is of type GCTAResponseTable. The
 * energy axis will be set to log10, the offset angle axis to radians.
 ***************************************************************************/
void GCTAAeff2D::read(const GFits* fits)
{
    // Clear instance
    clear();

    // Get pointer to effective area HDU
    GFitsTable* hdu = fits->table("EFFECTIVE AREA");

    // Read effective area table
    m_aeff.read(hdu);

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
 * @brief Write effective area into FITS file
 *
 * @param[in] fits FITS file.
 *
 * @todo Implement method.
 ***************************************************************************/
void GCTAAeff2D::write(GFits& fits) const
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Print effective area information
 *
 * @return Content of effective area instance.
 ***************************************************************************/
std::string GCTAAeff2D::print(void) const
{
    // Initialise result string
    std::string result;

    // Compute energy boundaries in TeV
    double emin = m_aeff.axis_lo(0,0);
    double emax = m_aeff.axis_hi(0,m_aeff.axis(0)-1);

    // Compute offset angle boundaries in deg
    double omin = m_aeff.axis_lo(1,0);
    double omax = m_aeff.axis_hi(1,m_aeff.axis(1)-1);

    // Append header
    result.append("=== GCTAAeff2D ===");
    result.append("\n"+parformat("Number of energy bins")+str(m_aeff.axis(0)));
    result.append("\n"+parformat("Number of offset bins")+str(m_aeff.axis(1)));
    result.append("\n"+parformat("Log10(Energy) range"));
    result.append(str(emin)+" - "+str(emax)+" TeV");
    result.append("\n"+parformat("Offset angle range"));
    result.append(str(omin)+" - "+str(omax)+" deg");

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
    m_aeff = aeff.m_aeff;

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
