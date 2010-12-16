/***************************************************************************
 *                 GLATAeff.cpp  -  Fermi LAT effective area               *
 * ----------------------------------------------------------------------- *
 *  copyright : (C) 2008-2010 by Jurgen Knodlseder                         *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GLATAeff.cpp
 * @brief Fermi LAT effective area class definition.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GLATAeff.hpp"
#include "GLATException.hpp"
#include "GTools.hpp"
#include "GFitsTableCol.hpp"
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
 * @param[in] logE Logarithm of the true photon energy (MeV).
 * @param[in] ctheta Cosine of the zenith angle with respect to the pointing 
 *            axis.
 ***************************************************************************/
double GLATAeff::operator() (const double& logE, const double& ctheta)
{
    // Get effective area value
    double aeff = (ctheta >= m_min_ctheta) 
                  ? m_aeff_bins.interpolate(logE, ctheta, m_aeff) : 0.0;

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
 * @brief Set minimum cos(theta) angle for effective area access
 *
 * @param[in] ctheta Cosine of the maximum zenith angle for which effective
 *                   areas will be returned (0 is returned for larger values)
 ***************************************************************************/
void GLATAeff::costhetamin(const double& ctheta)
{
    // Set minimum cos(theta) value
    m_min_ctheta = ctheta;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set livetime cube energy
 *
 * @param[in] energy Livetime cube energy.
 *
 * The livetime cube energy is used by the ltcube_ctheta() and
 * ltcube_ctheta_phi() methods that need an energy. Since we cannot pass
 * the energy as a method parameter we save it as a member of the class
 * so that it is avaible to the methods.
 ***************************************************************************/
void GLATAeff::ltcube_energy(const GEnergy& energy)
{
    // Set energy
    m_ltcube_logE = log10(energy.MeV());

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns the effective area as function of cos(theta) for the
 *        livetime cube energy
 *
 * @param[in] costheta Cosine of zenith angle.
 ***************************************************************************/
double GLATAeff::ltcube_ctheta(const double& costheta)
{
    // Return
    return (*this)(m_ltcube_logE, costheta);
}


/***********************************************************************//**
 * @brief Returns the effective area as function of cos(theta) for the
 *        livetime cube energy
 *
 * @param[in] costheta Cosine of zenith angle.
 * @param[in] phi Azimuth angle.
 *
 * @todo Not yet implemented.
 ***************************************************************************/
double GLATAeff::ltcube_ctheta_phi(const double& costheta, const double& phi)
{
    // Return
    return (*this)(m_ltcube_logE, costheta);
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
    m_aeff       = NULL;
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
    m_min_ctheta = aeff.m_min_ctheta;

    // Copy effective area array
    if (aeff.m_aeff != NULL) {
        int size = m_aeff_bins.size();
        if (size > 0) {
            m_aeff = new double[size];
            for (int i = 0; i < size; ++i)
                m_aeff[i] = aeff.m_aeff[i];
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GLATAeff::free_members(void)
{
    // Free effective area memory
    if (m_aeff != NULL) delete [] m_aeff;

    // Signal that effective area memory is free
    m_aeff = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read effective area from FITS table
 *
 * @param[in] hdu FITS table pointer.
 *
 * @exception GException::fits_column_not_found
 *            Effective area column not found
 * @exception GLATException::inconsistent_response
 *            Inconsistent response table encountered
 *
 * The effective area is converted into units of cm2.
 ***************************************************************************/
void GLATAeff::read_aeff(const GFitsTable* hdu)
{
    // Get energy and cos theta bins in response table
    m_aeff_bins.read(hdu);

    // Set minimum cos(theta)
    m_min_ctheta = m_aeff_bins.costheta_lo(0);

    // Continue only if there are effective area bins
    int size = m_aeff_bins.size();
    if (size > 0) {

        // Free effective area memory
        if (m_aeff != NULL) delete [] m_aeff;

        // Allocate memory for effective area
        m_aeff = new double[size];

        // Get pointer to effective area column
        GFitsTableCol* ptr = ((GFitsTable*)hdu)->column("EFFAREA");
        if (ptr == NULL)
            throw GException::fits_column_not_found(G_READ_AEFF, "EFFAREA");

        // Check consistency of effective area table
        int num = ptr->number();
        if (num != size)
            throw GLATException::inconsistent_response(G_READ_AEFF, num, size);

        // Copy data and convert from m2 into cm2
        for (int i = 0; i < size; ++i)
            m_aeff[i] = ptr->real(0,i) * 1.0e4;

    } // endif: there were effective area bins

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
