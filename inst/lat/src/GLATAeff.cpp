/***************************************************************************
 *                  GLATAeff.cpp - Fermi LAT effective area                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2016 by Juergen Knoedlseder                         *
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
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GException.hpp"
#include "GFilename.hpp"
#include "GEnergy.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"
#include "GFitsBinTable.hpp"
#include "GFitsTableCol.hpp"
#include "GFitsTableFloatCol.hpp"
#include "GLATAeff.hpp"
#include "GLATException.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_READ_AEFF                        "GLATAeff::read_aeff(GFitsTable&)"

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
 * @param[in] evtype Event type.
 *
 * Construct instance by loading the information from a FITS file. Both the
 * effective area information and the efficiency factor parameters are loaded
 * when available.
 ***************************************************************************/
GLATAeff::GLATAeff(const GFilename& filename, const std::string& evtype)
{
    // Initialise class members
    init_members();

    // Load effective area from file
    load(filename, evtype);

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
 * @return Effective area.
 ***************************************************************************/
GLATAeff& GLATAeff::operator=(const GLATAeff& aeff)
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
 *
 * Returns the effective area in units of cm2 for a given energy and
 * cos(theta) angle. The effective area is bi-linearily interpolated in the
 * log10(energy) - cos(theta) plane. The method assures that the effective
 * area value never becomes negative.
 ***************************************************************************/
double GLATAeff::operator()(const double& logE, const double& ctheta)
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
 * Returns the effective area in units of cm2 for a given energy, cos(theta)
 * angle and azimuth angle. The effective area is bi-linearily interpolated
 * in the log10(energy) - cos(theta) plane. The method assures that the
 * effective area value never becomes negative.
 *
 * @todo Phi-dependence not yet implemented.
 ***************************************************************************/
double GLATAeff::operator()(const double& logE, const double& ctheta,
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


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear effective area
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
 * @brief Clone effective area
 *
 * @return Pointer to deep copy of effective area.
 ***************************************************************************/
GLATAeff* GLATAeff::clone(void) const
{
    return new GLATAeff(*this);
}


/***********************************************************************//**
 * @brief Load effective area from FITS file
 *
 * @param[in] filename FITS file.
 * @param[in] evtype Event type.
 *
 * This method loads the effective area information, and if available, the
 * efficiency factors, from the FITS response file. See the GLATAeff::read
 * method for details.
 ***************************************************************************/
void GLATAeff::load(const GFilename& filename, const std::string& evtype)
{
    // Store event type
    m_evtype = evtype;

    // Open FITS file
    GFits fits(filename);

    // Read effective area from file
    read(fits);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save effective area into FITS file
 *
 * @param[in] filename FITS file.
 * @param[in] clobber Overwrite existing file? (default: false)
 *
 * This method saves the effective area information, and if available, the
 * efficiency factors, into the FITS response file. See the GLATAeff::write
 * method for details.
 ***************************************************************************/
void GLATAeff::save(const GFilename& filename, const bool& clobber)
{
    // Create FITS file
    GFits fits;

    // Write effective area into file
    write(fits);

    // Close FITS file
    fits.saveto(filename, clobber);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read effective area from FITS file
 *
 * @param[in] fits FITS file.
 *
 * Reads the effective area and efficiency parameter information form the
 * FITS file. The effective area is read from the extension EFFECTIVE AREA
 * (or EFFECTIVE AREA_<evtype> for Pass 8), the efficiency parameter
 * information from the extension EFFICIENCY_PARAMS (or 
 * EFFICIENCY_PARAMS_<evtype> for Pass 8). If the latter extension does not
 * exist, no efficiency parameters will be loaded.
 *
 * @todo Implement reading of Phi-dependence information.
 ***************************************************************************/
void GLATAeff::read(const GFits& fits)
{
    // Clear instance (keep event type)
    std::string evtype = m_evtype;
    clear();
    m_evtype = evtype;

    // Set extension names
    std::string effarea  = "EFFECTIVE AREA";
    std::string effparms = "EFFICIENCY_PARAMS";
    if (!fits.contains(effarea)) {
        effarea += "_" + m_evtype;
    }
    if (!fits.contains(effparms)) {
        effparms += "_" + m_evtype;
    }

    // Get pointer to effective area HDU
    const GFitsTable& hdu_aeff = *fits.table(effarea);

    // Read effective area
    read_aeff(hdu_aeff);

    // Read efficiency factors from appropriate HDU. Does nothing if the
    // HDU does not exist.
    if (fits.contains(effparms)) {
        const GFitsTable& hdu_eff = *fits.table(effparms);
        read_efficiency(hdu_eff);
    }

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

    // Write efficiency factors
    write_efficiency(fits);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Signals whether efficiency factors are present
 *
 * Verifies whether 2 efficiency factors are present. If efficiency factors
 * are included in the effective area FITS file, they will be automatically
 * loaded from the file, and efficiency factor functors will be allocated.
 * Otherwise, the functors point towards NULL.
 *
 * The method returns true if both efficiency factor functors exist, false
 * otherwise.
 ***************************************************************************/
bool GLATAeff::has_efficiency(void) const
{
    // Return
    return (m_eff_func1 != NULL && m_eff_func2 != NULL);
}


/***********************************************************************//**
 * @brief Returns efficiency factor 1
 *
 * @param[in] srcEng True energy of photon.
 *
 * Returns the efficiency factor 1 as function of the true photon energy.
 *
 * If no efficiency factors are present returns 1.
 *
 * @todo Implement cache to save computation time if called with same energy
 *       value (happens for binned analysis for example)
 ***************************************************************************/
double GLATAeff::efficiency_factor1(const GEnergy& srcEng) const
{
    // Initialise factor
    double factor = 1.0;

    // Compute efficiency factor. Note that the factor 1 uses the functor 2,
    // following the philosophy implemented in the ScienceTools method
    // EfficiencyFactor::getLivetimeFactors
    if (m_eff_func2 != NULL) {
        factor = (*m_eff_func2)(srcEng.log10MeV());
    }

    // Return factor
    return factor;
}


/***********************************************************************//**
 * @brief Returns efficiency factor 2
 *
 * @param[in] srcEng True energy of photon.
 *
 * Returns the efficiency factor 2 as function of the true photon energy.
 *
 * If no efficiency factors are present returns 2.
 *
 * @todo Implement cache to save computation time if called with same energy
 *       value (happens for binned analysis for example)
 ***************************************************************************/
double GLATAeff::efficiency_factor2(const GEnergy& srcEng) const
{
    // Initialise factor
    double factor = 0.0;

    // Compute efficiency factor. Note that the factor 2 uses the functor 1,
    // following the philosophy implemented in the ScienceTools method
    // EfficiencyFactor::getLivetimeFactors
    if (m_eff_func1 != NULL) {
        factor = (*m_eff_func1)(srcEng.log10MeV());
    }

    // Return factor
    return factor;
}


/***********************************************************************//**
 * @brief Print effective area information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing effective area information.
 ***************************************************************************/
std::string GLATAeff::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GLATAeff ===");

        // Append information
        result.append("\n"+gammalib::parformat("Number of energy bins") +
                      gammalib::str(nenergies()));
        result.append("\n"+gammalib::parformat("Number of cos theta bins") +
                      gammalib::str(ncostheta()));
        result.append("\n"+gammalib::parformat("Detector section")+m_evtype);
        result.append("\n"+gammalib::parformat("Efficiency factors"));
        if (has_efficiency()) {
            result.append("present");
        }
        else {
            result.append("absent");
        }

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
void GLATAeff::init_members(void)
{
    // Initialise members
    m_evtype.clear();
    m_aeff_bins.clear();
    m_aeff.clear();
    m_min_ctheta = 0.0;
    m_eff_func1  = NULL;
    m_eff_func2  = NULL;

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
    m_evtype     = aeff.m_evtype;
    m_aeff_bins  = aeff.m_aeff_bins;
    m_aeff       = aeff.m_aeff;
    m_min_ctheta = aeff.m_min_ctheta;

    // Clone functors
    m_eff_func1 = (aeff.m_eff_func1 != NULL) ? aeff.m_eff_func1->clone() : NULL;
    m_eff_func2 = (aeff.m_eff_func2 != NULL) ? aeff.m_eff_func2->clone() : NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GLATAeff::free_members(void)
{
    // Free functors
    if (m_eff_func1 != NULL) delete m_eff_func1;
    if (m_eff_func2 != NULL) delete m_eff_func2;

    // Signal that no functors are allocated
    m_eff_func1 = NULL;
    m_eff_func2 = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read effective area from FITS table
 *
 * @param[in] hdu FITS table.
 *
 * @exception GLATException::inconsistent_response
 *            Inconsistent response table encountered
 *
 * The effective area is converted into units of cm2.
 ***************************************************************************/
void GLATAeff::read_aeff(const GFitsTable& hdu)
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
        const GFitsTableCol* ptr = hdu["EFFAREA"];

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

    // Set detector section using the DETNAM keyword in the HDU
    /*
    std::string detnam = gammalib::strip_whitespace(hdu.string("DETNAM"));
    m_front            = (detnam == "FRONT");
    m_back             = (detnam == "BACK");
    */

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read efficiency factor parameters from FITS table
 *
 * @param[in] hdu FITS table.
 *
 * Reads the efficiency factor parameters from the EFFICIENCY_PARS column of
 * the EFFICIENCY_PARAMS extension in the effective area file. Note that the
 * column contains efficiency factor parameters for the front and back
 * section of the LAT detector, and we read here only those factors that
 * apply to the section for which the effective area has been defined. In
 * that way, we keep the efficiency factors separate for each detector
 * section, which is more precise than the handling in ScienceTools which
 * apply average efficiency factors.
 *
 * This method does not read efficiency factor parameters if neither the
 * front nor the back section of the detector has been identified.
 ***************************************************************************/
void GLATAeff::read_efficiency(const GFitsTable& hdu)
{
    // Free any existing functors
    if (m_eff_func1 != NULL) delete m_eff_func1;
    if (m_eff_func2 != NULL) delete m_eff_func2;

    // Signal that no functors are allocated
    m_eff_func1 = NULL;
    m_eff_func2 = NULL;

    // Get pointer to efficiency factors column
    const GFitsTableCol* ptr = hdu["EFFICIENCY_PARS"];

    // Allocate vectors to hold the parameters
    std::vector<double> par1;
    std::vector<double> par2;

    // If we have 2 rows in the table we have a Pass 8 response and hence
    // the table only applies to the specific event type
    if (ptr->length() == 2) {
        for (int i = 0; i < 6; ++i) {
            par1.push_back(ptr->real(0,i));
            par2.push_back(ptr->real(1,i));
        }
    }

    // ... otherwise we have an earlier response, and the front parameters
    // are stored in the first 2 rows while the back parameters are stored
    // in the second 2 rows
    else {
    
        // If we have a front section then read the front parameters
        if (m_evtype == "FRONT") {
            for (int i = 0; i < 6; ++i) {
                par1.push_back(ptr->real(0,i));
                par2.push_back(ptr->real(1,i));
            }
        }

        // If we have a back section then read the back parameters
        else if (m_evtype == "BACK") {
            for (int i = 0; i < 6; ++i) {
                par1.push_back(ptr->real(2,i));
                par2.push_back(ptr->real(3,i));
            }
        }
    }

    // If we have parameters, then allocate the efficiency factor
    // functors now
    if (par1.size() == 6) {
        m_eff_func1 = new GLATEfficiency(par1);
    }
    if (par2.size() == 6) {
        m_eff_func2 = new GLATEfficiency(par2);
    }

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
        m_aeff_bins.write(*hdu_aeff);

        // Allocate floating point vector columns
        GFitsTableFloatCol col_aeff = GFitsTableFloatCol("EFFAREA",  1, size);

        // Fill columns
        for (int i = 0; i < size; ++i) {
            col_aeff(0,i) = m_aeff[i];
        }

        // Append columns to table
        hdu_aeff->append(col_aeff);

        // Append HDU to FITS file
        file.append(*hdu_aeff);

        // Free binary table
        delete hdu_aeff;

    } // endif: there were data to write

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write efficiency factors into FITS file
 *
 * @param[in] file FITS file.
 *
 * @todo To be implemented.
 ***************************************************************************/
void GLATAeff::write_efficiency(GFits& file) const
{
    // Return
    return;
}
