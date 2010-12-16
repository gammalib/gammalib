/***************************************************************************
 *                GLATLtCube.cpp  -  Fermi LAT lifetime cube               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GLATLtCube.hpp"

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
GLATLtCube::GLATLtCube(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief File constructor
 *
 * @param[in] filename Lifetime cube filename.
 ***************************************************************************/
GLATLtCube::GLATLtCube(const std::string& filename)
{
    // Initialise class members
    init_members();

    // Load lifetime cube from file
    load(filename);

    // Return
    return;
}




/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] cube Lifetime cube.
 ***************************************************************************/
GLATLtCube::GLATLtCube(const GLATLtCube& cube)
{
    // Initialise class members
    init_members();

    // Copy members
    copy_members(cube);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GLATLtCube::~GLATLtCube(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Operators                               =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] cube Lifetime cube.
 ***************************************************************************/
GLATLtCube& GLATLtCube::operator= (const GLATLtCube& cube)
{
    // Execute only if object is not identical
    if (this != &cube) {

        // Free members
        free_members();

        // Initialise private members
        init_members();

        // Copy members
        copy_members(cube);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/***********************************************************************//**
 * @brief Sum function multiplied by efficiency corrected livetime over
 *        zenith angle
 *
 * @param[in] dir Sky direction.
 * @param[in] energy Energy.
 * @param[in] fct Function to evaluate.
 *
 * Computes
 * \f[\sum_{\cos \theta} T_{\rm corr.~live}(\cos \theta) f(\cos \theta)\f]
 * where
 * \f$T_{\rm corr.~live}(\cos \theta)\f$ is the efficiency corrected
 * livetime as a function of the cosine of the zenith angle, and
 * \f$f(\cos \theta)\f$ is a function that depends on the cosine of the
 * zenith angle.
 *
 * @todo Implement computation of efficiency factors
 ***************************************************************************/
double GLATLtCube::operator() (const GSkyDir& dir, const GEnergy& energy,
                               _ltcube_ctheta fct)
{
    // Compute exposure
    double exposure = m_exposure(dir, fct);

    // Optionally compute livetime factors for trigger rate- and
    // energy-dependent efficiency corrections
    if (m_livetime_correct) {
    
        // Compute correction factors
        double f1 = 1.0;
        double f2 = 0.0;
        //m_efficiencyFactor->getLivetimeFactors(energy, f1, f2);

        // Compute correction
        double correction = m_weighted_exposure(dir, fct);

        // Set exposure
        exposure = f1 * exposure + f2 * correction;

    } // endif: corrections requested

    // Return exposure
    return exposure;
}


/***********************************************************************//**
 * @brief Sum function multiplied by efficiency corrected livetime over
 *        zenith and azimuth angles
 *
 * @param[in] dir Sky direction.
 * @param[in] energy Energy.
 * @param[in] fct Function to evaluate.
 *
 * Computes
 * \f[\sum_{\cos \theta, \phi} T_{\rm corr.~live}(\cos \theta, \phi)
 *    f(\cos \theta, \phi)\f]
 * where
 * \f$T_{\rm corr.~live}(\cos \theta, \phi)\f$ is the efficiency corrected
 * livetime as a function of the cosine of the zenith and of the azimuth
 * angle, and
 * \f$f(\cos \theta, \phi)\f$ is a function that depends on the cosine of
 * the zenith angle and of the azimuth angle.
 *
 * @todo Implement computation of efficiency factors
 ***************************************************************************/
double GLATLtCube::operator() (const GSkyDir& dir, const GEnergy& energy,
                               _ltcube_ctheta_phi fct)
{
    // Compute exposure
    double exposure = m_exposure(dir, fct);

    // Optionally compute livetime factors for trigger rate- and
    // energy-dependent efficiency corrections
    if (m_livetime_correct) {
    
        // Compute correction factors
        double f1 = 1.0;
        double f2 = 0.0;
        //m_efficiencyFactor->getLivetimeFactors(energy, f1, f2);

        // Compute correction
        double correction = m_weighted_exposure(dir, fct);

        // Set exposure
        exposure = f1 * exposure + f2 * correction;

    } // endif: corrections requested

    // Return exposure
    return exposure;
}


/***********************************************************************//**
 * @brief Sum effective area multiplied by efficiency corrected livetime
 *        over zenith and (optionally) azimuth angles
 *
 * @param[in] dir Sky direction.
 * @param[in] energy Energy.
 * @param[in] aeff Effective area.
 *
 * Computes
 * \f[\sum_{\cos \theta, \phi} T_{\rm corr.~live}(\cos \theta, \phi)
 *    A_{\rm eff}(\cos \theta, \phi)\f]
 * where
 * \f$T_{\rm corr.~live}(\cos \theta, \phi)\f$ is the efficiency corrected
 * livetime as a function of the cosine of the zenith and of the azimuth
 * angle, and
 * \f$A_{\rm eff}(\cos \theta, \phi)\f$ is the effective area that depends
 * on the cosine of the zenith angle and (optionally) of the azimuth angle.
 *
 * @todo Implement computation of efficiency factors
 ***************************************************************************/
double GLATLtCube::operator() (const GSkyDir& dir, const GEnergy& energy,
                               const GLATAeff& aeff)
{
    // Compute exposure
    double exposure = m_exposure(dir, energy, aeff);

    // Optionally compute livetime factors for trigger rate- and
    // energy-dependent efficiency corrections
    if (m_livetime_correct) {
    
        // Compute correction factors
        double f1 = 1.0;
        double f2 = 0.0;
        //m_efficiencyFactor->getLivetimeFactors(energy, f1, f2);

        // Compute correction
        double correction = m_weighted_exposure(dir, energy, aeff);

        // Set exposure
        exposure = f1 * exposure + f2 * correction;

    } // endif: corrections requested

    // Return exposure
    return exposure;
}


/***********************************************************************//**
 * @brief Sum point spread function multiplied by efficiency corrected
 *        livetime over zenith angles
 *
 * @param[in] dir Sky direction.
 * @param[in] energy Energy.
 * @param[in] offset Offset from true direction (deg).
 * @param[in] psf Point spread function.
 *
 * Computes
 * \f[\sum_{\cos \theta, \phi} T_{\rm corr.~live}(\cos \theta, \phi)
 *    PSF(\log E, \delta, \cos \theta)\f]
 * where
 * \f$T_{\rm corr.~live}(\cos \theta, \phi)\f$ is the efficiency corrected
 * livetime as a function of the cosine of the zenith and of the azimuth
 * angle, and
 * \f$PSF(\log E, \delta, \cos \theta)\f$ is the point spread function that
 * depends on
 * the log10 of the energy (in MeV),
 * the offset angle from the true direction (in degrees), and
 * the cosine of the zenith angle.
 *
 * @todo Implement computation of efficiency factors
 ***************************************************************************/
double GLATLtCube::operator() (const GSkyDir& dir, const GEnergy& energy,
                               const double& offset, const GLATPsf& psf)
{
    // Compute exposure
    double exposure = m_exposure(dir, energy, offset, psf);

    // Optionally compute livetime factors for trigger rate- and
    // energy-dependent efficiency corrections
    if (m_livetime_correct) {
    
        // Compute correction factors
        double f1 = 1.0;
        double f2 = 0.0;
        //m_efficiencyFactor->getLivetimeFactors(energy, f1, f2);

        // Compute correction
        double correction = m_weighted_exposure(dir, energy, offset, psf);

        // Set exposure
        exposure = f1 * exposure + f2 * correction;

    } // endif: corrections requested

    // Return exposure
    return exposure;
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
void GLATLtCube::clear(void)
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
GLATLtCube* GLATLtCube::clone(void) const
{
    return new GLATLtCube(*this);
}


/***********************************************************************//**
 * @brief Load lifetime cube from FITS file
 *
 * @param[in] filename FITS file name.
 *
 * @todo Loading of cos theta boundaries not yet implemented.
 ***************************************************************************/
void GLATLtCube::load(const std::string& filename)
{
    // Clear object
    clear();

    // Open lifetime cube FITS file
    GFits file(filename);

    // Get HDUs
    GFitsTable* hdu_exposure          = file.table("EXPOSURE");
    GFitsTable* hdu_weighted_exposure = file.table("WEIGHTED_EXPOSURE");
    GFitsTable* hdu_cthetabounds      = file.table("CTHETABOUNDS");
    GFitsTable* hdu_gti               = file.table("GTI");

    // Load exposure
    m_exposure.read(hdu_exposure);

    // Load weighted exposure
    m_weighted_exposure.read(hdu_weighted_exposure);

    // Load cos theta boundaries
    
    // Load GTIs
    m_gti.read(hdu_gti);

    // Close FITS file
    file.close();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save lifetime cube into FITS file
 *
 * @param[in] filename FITS file name.
 * @param[in] clobber Overwrite existing file?
 *
 * @todo Not yet implemented.
 ***************************************************************************/
void GLATLtCube::save(const std::string& filename, bool clobber) const
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Print lifetime cube information
 ***************************************************************************/
std::string GLATLtCube::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GLATLtCube ===");
    result.append("\n"+m_exposure.print());
    result.append("\n"+m_weighted_exposure.print());
    //result.append("\n"+m_gti.print());

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GLATLtCube::init_members(void)
{
    // Initialise members
    m_exposure.clear();
    m_weighted_exposure.clear();
    m_gti.clear();
    m_livetime_correct = false;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] cube Lifetime cube.
 ***************************************************************************/
void GLATLtCube::copy_members(const GLATLtCube& cube)
{
    // Copy members
    m_exposure          = cube.m_exposure;
    m_weighted_exposure = cube.m_weighted_exposure;
    m_gti               = cube.m_gti;
    m_livetime_correct  = cube.m_livetime_correct;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GLATLtCube::free_members(void)
{
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
 * @param[in] cube Lifetime cube.
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GLATLtCube& cube)
{
     // Write lifetime cube in output stream
    os << cube.print();

    // Return output stream
    return os;
}


/***********************************************************************//**
 * @brief Log operator
 *
 * @param[in] log Logger.
 * @param[in] cube Lifetime cube.
 ***************************************************************************/
GLog& operator<< (GLog& log, const GLATLtCube& cube)
{
    // Write lifetime cube into logger
    log << cube.print();

    // Return logger
    return log;
}
