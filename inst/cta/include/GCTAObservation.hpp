/***************************************************************************
 *               GCTAObservation.hpp  -  CTA Observation class             *
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
/**
 * @file GCTAObservation.hpp
 * @brief GCTAObservation class interface definition.
 * @author J. Knodlseder
 */

#ifndef GCTAOBSERVATION_HPP
#define GCTAOBSERVATION_HPP

/* __ Includes ___________________________________________________________ */
#include "GObservation.hpp"


/***********************************************************************//**
 * @class GCTAObservation
 *
 * @brief Interface class for CTA observations.
 ***************************************************************************/
class GCTAObservation : public GObservation {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GCTAObservation& obs);

public:
    // Constructors and destructors
    GCTAObservation(void);
    GCTAObservation(const GCTAObservation& obs);
    virtual ~GCTAObservation(void);

    // Operators
    GCTAObservation& operator= (const GCTAObservation& obs);

    // Methods
    void response(const std::string& irfname, std::string caldb = "");
    void load_unbinned(const std::string& evname);
  
protected:
    // Protected methods
    void             init_members(void);
    void             copy_members(const GCTAObservation& obs);
    void             free_members(void);
    GCTAObservation* clone(void) const;

    // Protected data area
};

#endif /* GCTAOBSERVATION_HPP */
