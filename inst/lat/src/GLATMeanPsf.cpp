/***************************************************************************
 *                 GLATMeanPsf.cpp  -  Fermi LAT mean PSF                  *
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
#include <cmath>
#include "GLATMeanPsf.hpp"

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
GLATMeanPsf::GLATMeanPsf(void)
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] psf Mean PSF.
 ***************************************************************************/
GLATMeanPsf::GLATMeanPsf(const GLATMeanPsf& psf)
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
GLATMeanPsf::~GLATMeanPsf(void)
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
 * @param[in] psf Mean PSF.
 ***************************************************************************/
GLATMeanPsf& GLATMeanPsf::operator= (const GLATMeanPsf& psf)
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

/***********************************************************************//**
 * @brief Clear instance
 *
 * This method properly resets the object to an initial state.
 ***************************************************************************/
void GLATMeanPsf::clear(void)
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
GLATMeanPsf* GLATMeanPsf::clone(void) const
{
    return new GLATMeanPsf(*this);
}


/***********************************************************************//**
 * @brief Return size of mean PSF
***************************************************************************/
int GLATMeanPsf::size(void) const
{
    // Compute size
    int size = m_energy.size()*m_offset.size();

    // Return size
    return size;
}


/***********************************************************************//**
 * @brief Print lifetime cube information
 ***************************************************************************/
std::string GLATMeanPsf::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GLATMeanPsf ===");
    //result.append("\n"+m_dir.print());
    //result.append("\n"+m_weighted_exposure.print());
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
void GLATMeanPsf::init_members(void)
{
    // Initialise members
    m_dir.clear();
    m_psf.clear();
    m_exposure.clear();
    m_energy.clear();
    m_offset.clear();

    // Set offset array
    set_offsets();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] psf Mean PSF.
 ***************************************************************************/
void GLATMeanPsf::copy_members(const GLATMeanPsf& psf)
{
    // Copy members
    m_dir      = psf.m_dir;
    m_psf      = psf.m_psf;
    m_exposure = psf.m_exposure;
    m_energy   = psf.m_energy;
    m_offset   = psf.m_offset;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GLATMeanPsf::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set array of offset values in degrees
 *
 * The array of offset values defines the points at with the mean PSF will
 * be evaluated and stored.
 ***************************************************************************/
void GLATMeanPsf::set_offsets(void)
{
    // Set array parameters
    const double offset_min = 1.0e-4;
    const double offset_max = 70.0;
    const int    offset_num = 200;

    // Clear offset array
    m_offset.clear();
    m_offset.reserve(offset_num);

    // Set array
    double step = log(offset_max/offset_min)/(offset_num - 1.0);
    for (int i = 0; i < offset_num; ++i)
        m_offset.push_back(offset_min*exp(i*step));

    // Return
    return;
}


/***********************************************************************//**
 * @brief Compute mean PSF and exposure
 *
 * @param[in] obs LAT observation.
 ***************************************************************************/
void GLATMeanPsf::set_psf(const GLATObservation& obs)
{
    // Clear PSF and exposure arrays
    m_psf.clear();
    m_exposure.clear();

    // Allocate room for arrays
    m_psf.reserve(size());
    m_exposure.reserve(m_energy.size());

    // Loop over energies
    for (int ieng = 0; ieng < m_energy.size(); ++ieng) {

        // Compute exposure
        double exposure = 0.0;

        // Set exposure
        m_exposure.push_back(exposure);

        // Loop over all offsets
        for (int ioffset = 0; ioffset < m_offset.size(); ++ioffset) {

            // Initialise PSF
            double psf = 0.0;

            // Set PSF value
            m_psf.push_back(psf);

        } // endfor: looped over offsets
    } // endfor: looped over energies

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
 * @param[in] psf Mean PSF.
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GLATMeanPsf& psf)
{
     // Write mean PSF in output stream
    os << psf.print();

    // Return output stream
    return os;
}


/***********************************************************************//**
 * @brief Log operator
 *
 * @param[in] log Logger.
 * @param[in] psf Mean PSF.
 ***************************************************************************/
GLog& operator<< (GLog& log, const GLATMeanPsf& psf)
{
    // Write mean PSF into logger
    log << psf.print();

    // Return logger
    return log;
}
