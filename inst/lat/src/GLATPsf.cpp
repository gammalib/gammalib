/***************************************************************************
 *              GLATPsf.cpp - Fermi LAT point spread function              *
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
 * @file GLATPsf.cpp
 * @brief Fermi LAT point spread function class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GTools.hpp"
#include "GFilename.hpp"
#include "GFits.hpp"
#include "GFitsBinTable.hpp"
#include "GFitsTableFloatCol.hpp"
#include "GLATPsf.hpp"
#include "GLATPsfV1.hpp"
#include "GLATPsfV3.hpp"
#include "GLATException.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_READ                                        "GLATPsf::read(GFits&)"
#define G_WRITE                                      "GLATPsf::write(GFits&)"

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
 * Constructs an empty Fermi LAT point spread function.
 ***************************************************************************/
GLATPsf::GLATPsf(void)
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
 * Constructs a Fermi LAT point spread function by loading the point spread
 * function for a given event type from a FITS response file.
 ***************************************************************************/
GLATPsf::GLATPsf(const GFilename& filename, const std::string& evtype)
{
    // Initialise class members
    init_members();

    // Load PSF from FITS file
    load(filename, evtype);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] psf Point spread function.
 *
 * Constructs point spread function by copying point spread function from
 * another object.
 ***************************************************************************/
GLATPsf::GLATPsf(const GLATPsf& psf)
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
GLATPsf::~GLATPsf(void)
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
 * @return Point spread function.
 ***************************************************************************/
GLATPsf& GLATPsf::operator=(const GLATPsf& psf)
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


/***********************************************************************//**
 * @brief Return point spread function value
 *
 * @param[in] offset Offset angle (deg).
 * @param[in] logE Log10 of the true photon energy (MeV).
 * @param[in] ctheta Cosine of zenith angle.
 *
 * Returns the PSF value as function of the offset angle, the base 10
 * logarithm of the energy, and the cosine of the zenith angle. This method
 * calls the version dependent method.
 *
 * Returns 0 is no PSF has been allocated.
 ***************************************************************************/
double GLATPsf::operator()(const double& offset, const double& logE,
                           const double& ctheta)
{
    // Get PSF value
    double psf = (m_psf != NULL) ? m_psf->psf(offset, logE, ctheta) : 0.0;

    // Return PSF value
    return psf;
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
void GLATPsf::clear(void)
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
 *
 * @return Pointer to deep copy of point spread function.
 ***************************************************************************/
GLATPsf* GLATPsf::clone(void) const
{
    return new GLATPsf(*this);
}


/***********************************************************************//**
 * @brief Load point spread function from FITS file
 *
 * @param[in] filename FITS file.
 * @param[in] evtype Event type.
 *
 * Loads Fermi/LAT point spread function from FITS file.
 ***************************************************************************/
void GLATPsf::load(const GFilename& filename, const std::string& evtype)
{
    // Store event type
    m_evtype = evtype;

    // Open FITS file
    GFits fits(filename);

    // Read point spread function from file
    read(fits);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save point spread function into FITS file
 *
 * @param[in] filename FITS file.
 * @param[in] clobber Overwrite existing file?
 *
 * Saves Fermi/LAT point spread function into FITS file.
 ***************************************************************************/
void GLATPsf::save(const GFilename& filename, const bool& clobber)
{
    // Create FITS file
    GFits fits;

    // Write point spread function into file
    write(fits);

    // Close FITS file
    fits.saveto(filename, clobber);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read point spread function from FITS file
 *
 * @param[in] fits FITS file.
 *
 * @exception GException::invalid_response
 *            Invalid response type or unsupported response version found.
 *
 * Reads the Fermi LAT point spread function from FITS file. The method
 * determines the PSF version from the information found in the FITS file
 * and allocates the proper PSF version class. It reads the PSF information
 * from the FITS file. 
 *
 * @todo Implement PSF version 2.
 ***************************************************************************/
void GLATPsf::read(const GFits& fits)
{
    // Clear instance (keep event type)
    std::string evtype = m_evtype;
    clear();
    m_evtype = evtype;

    // Set extension names
    std::string rpsf    = "RPSF";
    std::string scaling = "PSF_SCALING_PARAMS";
    if (!fits.contains(rpsf)) {
        rpsf += "_" + m_evtype;
    }
    if (!fits.contains(scaling)) {
        scaling += "_" + m_evtype;
    }

    // Get PSF parameters table
    const GFitsTable& hdu_rpsf = *fits.table(rpsf);

    // Get PSF scaling parameters table
    const GFitsTable& hdu_scale = *fits.table(scaling);

    // Determine PSF version (default version is version 1)
    int version = (hdu_rpsf.has_card("PSFVER")) ? hdu_rpsf.integer("PSFVER") : 1;

    // Determine PSF type
    bool        front  = true;
    std::string detnam = 
        gammalib::strip_whitespace(gammalib::toupper(hdu_rpsf.string("DETNAM")));
    if (detnam == "FRONT") {
        front = true;
    }
    else if (detnam == "BACK") {
        front = false;
    }
    else {
        front = false; // If we are here we should have a Pass 8 response.
                       // For Pass 8 the setting of this flag is irrelevant.
    }

    // Allocate point spread function
    switch (version) {
    case 1:
        m_psf = new GLATPsfV1;
        break;
    case 3:
        m_psf = new GLATPsfV3;
        break;
    default:
        throw GLATException::invalid_response(G_READ, 
              "Unsupported response function version "+gammalib::str(version)+".");
        break;
    }

    // Set PSF attributes (needed before reading scaling parameters for
    // Pass 6 and Pass 7 response functions)
    m_psf->front(front);

    // Read PSF scaling parameters
    m_psf->read_scale(hdu_scale);

    // Read PSF
    m_psf->read(hdu_rpsf);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write point spread function into FITS file
 *
 * @param[in] fits FITS file.
 *
 * Writes the version dependent point spread function into the FITS file. The
 * method does nothing if no PSF has been allocated.
 *
 * @todo Implement PSF versions 2 and 3.
 ***************************************************************************/
void GLATPsf::write(GFits& fits) const
{
    // Continue only if PSF is valid
    if (m_psf != NULL) {

        // Write point spread function
        m_psf->write(fits);

        // Write scaling parameters
        m_psf->write_scale(fits);

    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns size of PSF
 *
 * @return Size of point spread function.
 *
 * Returns the size of the PSF which is defined as the number of energy bins
 * times the number of cos(theta) bins. If no PSF has been allocated, 0 is
 * returned.
 ***************************************************************************/
int GLATPsf::size(void) const
{
    // Return size of PSF
    return (nenergies()*ncostheta());
}    


/***********************************************************************//**
 * @brief Returns number of energy bins in PSF
 *
 * @return Number of energy bins.
 *
 * Returns the number of energy bins in PSF. If no PSF has been allocated,
 * 0 is returned.
 ***************************************************************************/
int GLATPsf::nenergies(void) const
{
    // Retrieve number of energy bins
    int nenergies = (m_psf != NULL) ? m_psf->nenergies() : 0;
    
    // Return number of energy bins
    return nenergies;
}    


/***********************************************************************//**
 * @brief Returns number of cos(theta) bins in PSF
 *
 * @return Number of cos(theta) bins.
 *
 * Returns the number of cos(theta) bins in PSF. If no PSF has been
 * allocated, 0 is returned.
 ***************************************************************************/
int GLATPsf::ncostheta(void) const
{
    // Retrieve number of cos(theta) bins
    int ncostheta = (m_psf != NULL) ? m_psf->ncostheta() : 0;
    
    // Return number of cos(theta) bins
    return ncostheta;
}    
    

/***********************************************************************//**
 * @brief Returns minimum cos(theta) angle
 *
 * @return Minimum cos(theta) angle.
 *
 * Returns the minimum cos(theta) angle for point spread function access. If
 * no PSF has been allocated, 0 is returned.
 ***************************************************************************/
double GLATPsf::costhetamin(void) const
{
    // Retrieve number of energy bins
    double costhetamin = (m_psf != NULL) ? m_psf->costhetamin() : 0.0;
    
    // Return minimum cos(theta)
    return costhetamin;
}    
    

/***********************************************************************//**
 * @brief Set minimum cos(theta) angle for point spread function access
 *
 * @param[in] ctheta Cosine of maximum zenith angle.
 *
 * Sets the minimum cos(theta) angle. No verification of the specified 
 * maximum cos(theta) value is done.
 ***************************************************************************/
void GLATPsf::costhetamin(const double& ctheta)
{
    // Continue only if PSF is valid
    if (m_psf != NULL) {

        // Set minimum cos(theta) value
        m_psf->costhetamin(ctheta);

    }
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns PSF version
 *
 * Returns the PSF version number.
 *
 * Returns 0 if no PSF has been allocated.
 ***************************************************************************/
int GLATPsf::version(void) const
{
    // Retrieve version number
    int version = (m_psf != NULL) ? m_psf->version() : 0;
    
    // Return version number
    return version;
}    


/***********************************************************************//**
 * @brief Print point spread function information
 *
 * @param[in] chatter Chattiness.
 * @return String containing point spread function information.
 ***************************************************************************/
std::string GLATPsf::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GLATPsf ===");
    
        // No PSF has been loaded ...
        if (m_psf == NULL) {
            result.append("\n"+gammalib::parformat("Version"));
            result.append("No PSF loaded");
        }
    
        // ... PSF has been loaded
        else {
            result.append("\n"+gammalib::parformat("Version"));
            result.append(gammalib::str(version()));
            result.append("\n"+gammalib::parformat("Event type"));
            result.append(evtype());
            result.append("\n"+gammalib::parformat("Energy scaling"));
            result.append("sqrt(");
            result.append("("+gammalib::str(m_psf->scale_par1())+"*(E/100)^");
            result.append(gammalib::str(m_psf->scale_index())+")^2");
            result.append(" + ("+gammalib::str(m_psf->scale_par2())+")^2)");
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
void GLATPsf::init_members(void)
{
    // Initialise members
    m_evtype.clear();
    m_psf = NULL;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] psf Point spread function.
 ***************************************************************************/
void GLATPsf::copy_members(const GLATPsf& psf)
{
    // Copy members
    m_evtype = psf.m_evtype;

    // Clone members
    if (psf.m_psf != NULL) {
        m_psf = psf.m_psf->clone();
    }
    else {
        m_psf = NULL;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GLATPsf::free_members(void)
{
    // Free PSF
    if (m_psf != NULL) delete m_psf;

    // Signal that PSF is free
    m_psf = NULL;

    // Return
    return;
}
