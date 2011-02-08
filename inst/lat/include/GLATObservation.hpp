/***************************************************************************
 *               GLATObservation.hpp  -  LAT Observation class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2011 by Jurgen Knodlseder                           *
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
 * @brief LAT Observation class interface definition
 * @author J. Knodlseder
 */

#ifndef GLATOBSERVATION_HPP
#define GLATOBSERVATION_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GObservation.hpp"
#include "GLATResponse.hpp"
#include "GLATPointing.hpp"
#include "GLATLtCube.hpp"
#include "GTime.hpp"
#include "GModel.hpp"
#include "GModels.hpp"


/***********************************************************************//**
 * @class GLATObservation
 *
 * @brief Interface for the LAT observation class
 ***************************************************************************/
class GLATObservation : public GObservation {

public:
    // Constructors and destructors
    GLATObservation(void);
    GLATObservation(const GLATObservation& obs);
    virtual ~GLATObservation(void);

    // Operators
    GLATObservation& operator= (const GLATObservation& obs);

    // Implemented pure virtual base class methods
    virtual void             clear(void);
    virtual GLATObservation* clone(void) const;
    virtual GLATResponse*    response(void) const;
    virtual GLATPointing*    pointing(const GTime& time) const;
    virtual std::string      instrument(void) const;
    virtual std::string      print(void) const;

    // Other methods
    void                     load_unbinned(const std::string& ft1name,
                                           const std::string& ft2name,
                                           const std::string& ltcube_name);
    void                     load_binned(const std::string& cntmap_name,
                                         const std::string& expmap_name,
                                         const std::string& ltcube_name);
    void                     response(const std::string& irfname,
                                      std::string caldb = "");    
    GLATLtCube*              ltcube(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GLATObservation& obs);
    void free_members(void);

    // Protected members
    GLATResponse* m_response;     //!< Pointer to instrument response functions
    GLATPointing* m_pointing;     //!< Pointer to pointing direction
    GLATLtCube*   m_ltcube;       //!< Pointer to lifetime cube
};

#endif /* GLATOBSERVATION_HPP */
