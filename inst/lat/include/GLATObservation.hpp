/***************************************************************************
 *               GLATObservation.hpp  -  LAT Observation class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2010 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GLATObservation.hpp
 * @brief GLATObservation class interface definition.
 * @author J. Knodlseder
 */

#ifndef GLATOBSERVATION_HPP
#define GLATOBSERVATION_HPP

/* __ Includes ___________________________________________________________ */
#include "GObservation.hpp"
#include "GLATResponse.hpp"
#include "GLATPointing.hpp"
#include "GTime.hpp"
#include "GModel.hpp"


/***********************************************************************//**
 * @class GLATObservation
 *
 * @brief Interface for the LAT observation classes.
 ***************************************************************************/
class GLATObservation : public GObservation {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GLATObservation& obs);

public:
    // Constructors and destructors
    GLATObservation();
    GLATObservation(const GLATObservation& obs);
    virtual ~GLATObservation();

    // Operators
    GLATObservation& operator= (const GLATObservation& obs);

    // Implement pure virtual methods
    GLATObservation* clone(void) const;
    void             response(const std::string& irfname, std::string caldb = "");
    GResponse*       response(const GTime& time) const;
    GPointing*       pointing(const GTime& time) const;
    std::string      instrument(void) const;

    // Other methods
    void load_unbinned(const std::string& ft1name, const std::string& ft2name,
                       const std::string& ltcube_name);
    void load_binned(const std::string& cntmap_name, const std::string& expmap_name,
                     const std::string& ltcube_name);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GLATObservation& obs);
    void free_members(void);

    // Protected data area
    GLATResponse* m_response;     //!< Pointer to instrument response functions
    GLATPointing* m_pointing;     //!< Pointer to pointing direction
};

#endif /* GLATOBSERVATION_HPP */
