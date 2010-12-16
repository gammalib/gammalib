/***************************************************************************
 *                  GLATMeanPsf.hpp  -  Fermi LAT mean PSF                 *
 * ----------------------------------------------------------------------- *
 *  copyright : (C) 2010 by Jurgen Knodlseder                              *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GLATMeanPsf.hpp
 * @brief GLATMeanPsf class definition.
 * @author J. Knodlseder
 */

#ifndef GLATMEANPSF_HPP
#define GLATMEANPSF_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include <iostream>
#include "GLog.hpp"
#include "GSkyDir.hpp"
#include "GEnergy.hpp"
#include "GLATObservation.hpp"


/***********************************************************************//**
 * @class GLATMeanPsf
 *
 * @brief Interface for the Fermi LAT position-dependent mean PSF.
 *
 * The position-dependent mean PSF is the point spread function that has
 * been averaged over the zenith and azimuth angles of an observation. The
 * averaging is done using the livetime cube which holds the lifetime as
 * function and zenith and azimuth angles for an observation.
 ***************************************************************************/
class GLATMeanPsf {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GLATMeanPsf& psf);
    friend GLog&         operator<< (GLog& log, const GLATMeanPsf& psf);

public:
    // Constructors and destructors
    GLATMeanPsf(void);
    GLATMeanPsf(const GSkyDir& dir, const GLATObservation& obs);
    GLATMeanPsf(const GLATMeanPsf& cube);
    virtual ~GLATMeanPsf(void);

    // Operators
    GLATMeanPsf& operator= (const GLATMeanPsf& cube);

    // Methods
    void         clear(void);
    GLATMeanPsf* clone(void) const;
    int          size(void) const;
    std::string  print(void) const;

private:
    // Methods
    void init_members(void);
    void copy_members(const GLATMeanPsf& psf);
    void free_members(void);
    void set_offsets(void);
    void set_psf(const GSkyDir& dir, const GLATObservation& obs);
    
    // Protected members
    GSkyDir              m_dir;        //!< Source direction for mean PSF
    std::vector<double>  m_psf;        //!< Mean PSF values
    std::vector<double>  m_exposure;   //!< Mean exposure
    std::vector<GEnergy> m_energy;     //!< Energies of mean PSF
    std::vector<double>  m_offset;     //!< Offsets of mean PSF
};

#endif /* GLATMEANPSF_HPP */
