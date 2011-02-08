/***************************************************************************
 *        GMWLObservation.i  -  Multi-wavelength observation class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GMWLObservation.i
 * @brief Multi-wavelength observation class Python interface definition
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GMWLObservation.hpp"
%}


/***********************************************************************//**
 * @class GMWLObservation
 *
 * @brief Python interface class for multi-wavelength observations
 *
 * This class implements a multi-wavelength observation. A multi-wavelength
 * observation contains spectral points obtained with an unspecified
 * instrument. The spectral points are given in physical units.
 ***************************************************************************/
class GMWLObservation : public GObservation {
public:
    // Constructors and destructors
    GMWLObservation(void);
    explicit GMWLObservation(const std::string& filename);
    explicit GMWLObservation(const std::string& filename, const std::string& extname);
    GMWLObservation(const GMWLObservation& obs);
    virtual ~GMWLObservation(void);

    // Implement pure virtual methods
    virtual void             clear(void);
    virtual GMWLObservation* clone(void) const;
    virtual GMWLResponse*    response(void) const;
    virtual GMWLPointing*    pointing(const GTime& time) const;
    virtual std::string      instrument(void) const;

    // Other methods
    void load(const std::string& filename);
    void load(const std::string& filename, const std::string& extname);
    void instrument(const std::string& instrument);
};


/***********************************************************************//**
 * @brief GMWLObservation class extension
 ***************************************************************************/
%extend GMWLObservation {
    GMWLObservation copy() {
        return (*self);
    }
};


/***********************************************************************//**
 * @brief GMWLObservation type casts
 ***************************************************************************/
%inline %{
    GMWLObservation* cast_GMWLObservation(GObservation* obs) {
        GMWLObservation* mwl = dynamic_cast<GMWLObservation*>(obs);
        if (mwl == NULL)
            throw GException::bad_type("cast_GMWLObservation(GObservation* obs)",
                                       "GObservation not of type GMWLObservation");            
        return mwl;
    }
%}
