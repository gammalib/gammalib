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
#include "GFits.hpp"


/***********************************************************************//**
 * @class GLATObservation
 *
 * @brief Interface for the LAT observation classes.
 ***************************************************************************/
class GLATObservation : public GObservation {

public:
    // Constructors and destructors
    GLATObservation();
    GLATObservation(const std::string& ft1name, const std::string& ft2name);
    GLATObservation(const GLATObservation& obs);
    virtual ~GLATObservation();

    // Operators
    GLATObservation& operator= (const GLATObservation& obs);

    // Methods
    GFits* ft1(void) const;
    GFits* ft2(void) const;
  
protected:
    // Protected methods
    void             init_members(void);
    void             copy_members(const GLATObservation& obs);
    void             free_members(void);
    GLATObservation* clone(void) const;

    // Protected data area

private:
    GFits* m_ft1;     //!< FT1 FITS file
    GFits* m_ft2;     //!< FT2 FITS file
    GFits* m_ltcube;  //!< Lifetime cube FITS file
};

#endif /* GLATOBSERVATION_HPP */
