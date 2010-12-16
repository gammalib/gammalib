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
#include "GLATAeff.hpp"
#include "GLATPsf.hpp"
#include "GLATException.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_SET                  "GLATMeanPsf::set(GSkyDir&, GLATObservation&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */

/* __ Constants __________________________________________________________ */

/* __ Typedefs ___________________________________________________________ */


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
 * @brief PSF constructor
 *
 * @param[in] dir Sky direction.
 * @param[in] obs LAT observation.
 ***************************************************************************/
GLATMeanPsf::GLATMeanPsf(const GSkyDir& dir, const GLATObservation& obs)
{
    // Initialise class members
    init_members();

    // Setup PSF
    set(dir, obs);

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
 * @brief Compute mean PSF and exposure
 *
 * @param[in] dir Source location.
 * @param[in] obs LAT observation.
 *
 * @exception GLATException::no_response
 *            Response has not been defined.
 * @exception GLATException::no_ltcube
 *            Livetime cube has not been defined.
 * @exception GLATException::no_ebds
 *            Energy binning has not been defined.
 *
 * Computes the mean PSF and the energy dependent exposure for a source at
 * a given sky location. The PSF is computed for all bin boundaries of the
 * observation, hence if there are N energy bins there will be N+1 energies
 * at which the mean PSF is computed.
 ***************************************************************************/
void GLATMeanPsf::set(const GSkyDir& dir, const GLATObservation& obs)
{
    // Clear PSF, exposure and energy arrays
    m_psf.clear();
    m_exposure.clear();
    m_energy.clear();

    // Get pointers on response, livetime cube and energy boundaries
    GLATResponse* rsp = obs.response();
    if (rsp == NULL)
        throw GLATException::no_response(G_SET);

    // Get pointer on livetime cube
    GLATLtCube* ltcube = obs.ltcube();
    if (ltcube == NULL)
        throw GLATException::no_ltcube(G_SET);
    
    // Get pointer on energy boundaries
    GEbounds* ebds = ((GLATObservation*)&obs)->ebounds();
    if (ebds == NULL)
        throw GLATException::no_ebds(G_SET);

    // Allocate room for arrays
    m_psf.reserve(size());
    m_exposure.reserve(m_energy.size());

    // Set energy nodes from the bin boundaries of the observations energy
    // boundaries.
    m_energy.reserve(ebds->size()+1);
    m_energy.push_back(ebds->emin(0));
    for (int i = 0; i < ebds->size(); ++i)
        m_energy.push_back(ebds->emax(i));

    // Loop over energies
    for (int ieng = 0; ieng < m_energy.size(); ++ieng) {

        // Compute exposure by looping over the responses
        double exposure = 0.0;
        for (int i = 0; i < rsp->size(); ++i)
            exposure += (*ltcube)(dir, m_energy[ieng], *rsp->aeff(i));

        // Set exposure
        m_exposure.push_back(exposure);

        // Loop over all offset angles
        for (int ioffset = 0; ioffset < m_offset.size(); ++ioffset) {

            // Compute point spread function by looping over the responses
            double psf = 0.0;
            for (int i = 0; i < rsp->size(); ++i)
                psf += (*ltcube)(dir, m_energy[ieng], m_offset[ioffset],
                                 *rsp->psf(i));

            // Normalize PSF by exposure and clip when exposure drops to 0
            psf = (exposure > 0) ? psf/exposure : 0.0;

            // Set PSF value
            m_psf.push_back(psf);

        } // endfor: looped over offsets
    } // endfor: looped over energies

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print lifetime cube information
 ***************************************************************************/
std::string GLATMeanPsf::print(void) const
{
    // Initialise result string
    std::string result;

    // Compute exposure range
    double min_exposure = 0.0;
    double max_exposure = 0.0;
    for (int i = 0; i < nenergies(); ++i) {
        if (i == 0) {
            min_exposure = m_exposure[i];
            max_exposure = m_exposure[i];
        }
        else {
            if (m_exposure[i] < min_exposure) min_exposure = m_exposure[i];
            if (m_exposure[i] > max_exposure) max_exposure = m_exposure[i];
        }
    }
    
    // Append header
    result.append("=== GLATMeanPsf ===");
    result.append("\n"+parformat("Source direction"));
    result.append("\n"+parformat("Offset angles")+str(noffsets()));
    result.append("\n"+parformat("Energy values")+str(nenergies()));
    result.append("\n"+parformat("Exposure range"));
    result.append(str(min_exposure)+" - "+str(max_exposure)+" s cm2");

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
