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


/***************************************************************************
 *                          GResponse class definition                     *
 ***************************************************************************/
/**
 * @class GResponse
 *
 * @brief Abstract interface for the instrument response function classes.
 *
 * @author J. Knodlseder
 */
class GResponse {

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
    /// @param obsDir Observed photon direction
    /// @param obsEng Observed energy of photon
    /// @param srcDir True photon direction
    /// @param srcEng True energy of photon
    /// @param instPntDir Instrument pointing direction (e.g. z-axis)
    /// @param instPosAng Instrument position angle
    /// @param time Photon arrival time
    virtual double irf(const GSkyDir& obsDir, const double& obsEng,
                       const GSkyDir& srcDir, const double& srcEng,
                       const GSkyDir& instPntDir, const double& instPosAng,
                       const double& time) = 0;

    /// Pure virtual method to define the interface for the member function
    /// returning the point spread function value
    /// @param obsDir Observed photon direction
    /// @param obsEng Observed energy of photon
    /// @param srcDir True photon direction
    /// @param srcEng True energy of photon
    /// @param instPntDir Instrument pointing direction (e.g. z-axis)
    /// @param instPosAng Instrument position angle
    /// @param time Photon arrival time
    virtual double psf(const GSkyDir& obsDir, const double& obsEng,
                       const GSkyDir& srcDir, const double& srcEng,
                       const GSkyDir& instPntDir, const double& instPosAng,
                       const double& time) = 0;

    /// Pure virtual method to define the interface for the member function
    /// returning the effective area value
    /// @param obsDir Observed photon direction
    /// @param obsEng Observed energy of photon
    /// @param srcDir True photon direction
    /// @param srcEng True energy of photon
    /// @param instPntDir Instrument pointing direction (e.g. z-axis)
    /// @param instPosAng Instrument position angle
    /// @param time Photon arrival time
    virtual double aeff(const GSkyDir& obsDir, const double& obsEng,
                        const GSkyDir& srcDir, const double& srcEng,
                        const GSkyDir& instPntDir, const double& instPosAng,
                        const double& time) = 0;

    /// Pure virtual method to define the interface for the member function
    /// returning the energy dispersion value
    /// @param obsDir Observed photon direction
    /// @param obsEng Observed energy of photon
    /// @param srcDir True photon direction
    /// @param srcEng True energy of photon
    /// @param instPntDir Instrument pointing direction (e.g. z-axis)
    /// @param instPosAng Instrument position angle
    /// @param time Photon arrival time
    virtual double edisp(const GSkyDir& obsDir, const double& obsEng,
                         const GSkyDir& srcDir, const double& srcEng,
                         const GSkyDir& instPntDir, const double& instPosAng,
                         const double& time) = 0;

    /// Pure virtual method to define the calibration database
    /// @param caldb Calibration database to be used
    virtual void set_caldb(std::string caldb) = 0;

    /// Pure virtual method to load the instrument response from the calibration
    /// database
    /// @param rspname Name of the instrument response
    virtual void load(std::string rspname) = 0;

    /// Pure virtual method to save the instrument response into the calibration
    /// database
    /// @param rspname Name of the instrument response
    virtual void save(std::string rspname) const = 0;

    /// Pure virtual method to clone the response instance
    virtual GResponse* clone(void) const = 0;

protected:
    /// Name of or path to the calibration database
    std::string m_caldb;

    /// Name of the instrument response
    std::string m_rspname;

    /// Protected method that initialises protected members of the class
    void init_members(void);

    /// Protected method that copies protected members of the class
    void copy_members(const GResponse& rsp);

    /// Protected method that frees memory allocated by instances of the class
    void free_members(void);

private:
};


/***************************************************************************
 *                              Inline methods                             *
 ***************************************************************************/

#endif /* GRESPONSE_HPP */
