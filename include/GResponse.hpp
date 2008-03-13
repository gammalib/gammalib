/***************************************************************************
 *              GResponse.hpp  -  Response abstract base class             *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2008 by Jurgen Knodlseder                   *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GResponse.hpp
 * @brief GResponse class definition.
 * @author J. Knodlseder
 */

#ifndef GRESPONSE_HPP
#define GRESPONSE_HPP

/* __ Includes ___________________________________________________________ */
#include "GSkyDir.hpp"

/* __ Namespaces _________________________________________________________ */


/***********************************************************************//**
 * @class GResponse
 *
 * @brief Abstract interface for the instrument response function classes.
 ***************************************************************************/
class GResponse {

  // Friend classes
  friend class GObservation;

public:
    /// Constructor
    GResponse();

    /// Copy constructor
    GResponse(const GResponse& rsp);

    /// Destructor
    virtual ~GResponse();

    /// Assignment operator
    virtual GResponse& operator= (const GResponse& rsp);

    /// Pure virtual method to define the interface for the member function
    /// returning the instrument response function value
    /// @param[in] obsDir Observed photon direction
    /// @param[in] obsEng Observed energy of photon
    /// @param[in] srcDir True photon direction
    /// @param[in] srcEng True energy of photon
    /// @param[in] instPntDir Instrument pointing direction (e.g. z-axis)
    /// @param[in] instPosAng Instrument position angle
    /// @param[in] time Photon arrival time
    virtual double irf(const GSkyDir& obsDir, const double& obsEng,
                       const GSkyDir& srcDir, const double& srcEng,
                       const GSkyDir& instPntDir, const double& instPosAng,
                       const double& time) = 0;

    /// Pure virtual method to define the interface for the member function
    /// returning the point spread function value
    /// @param[in] obsDir Observed photon direction
    /// @param[in] obsEng Observed energy of photon
    /// @param[in] srcDir True photon direction
    /// @param[in] srcEng True energy of photon
    /// @param[in] instPntDir Instrument pointing direction (e.g. z-axis)
    /// @param[in] instPosAng Instrument position angle
    /// @param[in] time Photon arrival time
    virtual double psf(const GSkyDir& obsDir, const double& obsEng,
                       const GSkyDir& srcDir, const double& srcEng,
                       const GSkyDir& instPntDir, const double& instPosAng,
                       const double& time) = 0;

    /// Pure virtual method to define the interface for the member function
    /// returning the effective area value
    /// @param[in] obsDir Observed photon direction
    /// @param[in] obsEng Observed energy of photon
    /// @param[in] srcDir True photon direction
    /// @param[in] srcEng True energy of photon
    /// @param[in] instPntDir Instrument pointing direction (e.g. z-axis)
    /// @param[in] instPosAng Instrument position angle
    /// @param[in] time Photon arrival time
    virtual double aeff(const GSkyDir& obsDir, const double& obsEng,
                        const GSkyDir& srcDir, const double& srcEng,
                        const GSkyDir& instPntDir, const double& instPosAng,
                        const double& time) = 0;

    /// Pure virtual method to define the interface for the member function
    /// returning the energy dispersion value
    /// @param[in] obsDir Observed photon direction
    /// @param[in] obsEng Observed energy of photon
    /// @param[in] srcDir True photon direction
    /// @param[in] srcEng True energy of photon
    /// @param[in] instPntDir Instrument pointing direction (e.g. z-axis)
    /// @param[in] instPosAng Instrument position angle
    /// @param[in] time Photon arrival time
    virtual double edisp(const GSkyDir& obsDir, const double& obsEng,
                         const GSkyDir& srcDir, const double& srcEng,
                         const GSkyDir& instPntDir, const double& instPosAng,
                         const double& time) = 0;

    /// Pure virtual method to define the calibration database
    /// @param[in] caldb Calibration database to be used
    virtual void set_caldb(const std::string& caldb) = 0;

protected:
    // Protected methods
    void    init_members(void);
    void    copy_members(const GResponse& rsp);
    void    free_members(void);
    virtual GResponse* clone(void) const = 0;
    
    // Protected data area 
    std::string m_caldb;    //!< Name of or path to the calibration database
    std::string m_rspname;  //!< Name of the instrument response



private:
};


/***************************************************************************
 *                              Inline methods                             *
 ***************************************************************************/

#endif /* GRESPONSE_HPP */
