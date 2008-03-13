/***************************************************************************
 *               GLATObservation.hpp  -  LAT Observation class             *
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
 * @file GLATObservation.hpp
 * @brief GLATObservation class interface definition.
 * @author J. Knodlseder
 */

#ifndef GLATOBSERVATION_HPP
#define GLATOBSERVATION_HPP

/* __ Includes ___________________________________________________________ */
#include "GObservation.hpp"

/* __ Namespaces _________________________________________________________ */


/***********************************************************************//**
 * @class GLATObservation
 *
 * @brief Interface for the LAT observation classes.
 ***************************************************************************/
class GLATObservation : public GObservation {

public:
    // Constructors and destructors
    GLATObservation();
    GLATObservation(const GLATObservation& obs);
    virtual ~GLATObservation();

    // Operators
    virtual GLATObservation& operator= (const GLATObservation& obs);

    // Methods
  
protected:
    // Protected methods
    void    init_members(void);
    void    copy_members(const GLATObservation& obs);
    void    free_members(void);
    virtual GLATObservation* clone(void) const;

    // Protected data area

private:
};

#endif /* GLATOBSERVATION_HPP */
