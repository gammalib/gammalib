/***************************************************************************
 *   GMWLObservation.i  -  Multi-wavelength observation class python I/F   *
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
 * @file GMWLObservation.i
 * @brief GMWLObservation class python interface
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GMWLObservation.hpp"
%}


/***********************************************************************//**
 * @class GMWLObservation
 *
 * @brief Interface class for multi-wavelength observations
 *
 * This class implements a multi-wavelength observation. A multi-wavelength
 * observation contains spectral points obtained with an unspecified
 * instrument. The spectral points are given in physical units.
 ***************************************************************************/
class GMWLObservation : public GObservation {
public:
    // Constructors and destructors
    GMWLObservation(void);
    GMWLObservation(const std::string& filename);
    GMWLObservation(const std::string& filename, const std::string& extname);
    GMWLObservation(const GMWLObservation& obs);
    virtual ~GMWLObservation(void);

    // Implement pure virtual methods
    GMWLObservation* clone(void) const;
    void             response(const std::string& irfname, std::string caldb = "");
    GResponse*       response(const GTime& time) const { return NULL; }
    GPointing*       pointing(const GTime& time) const { return NULL; }
    std::string      instrument(void) const { return m_instrument; }

    // Other methods
    void        load(const std::string& filename);
    void        load(const std::string& filename, const std::string& extname);
    void        instrument(const std::string& instrument) { m_instrument=instrument; }
};


/***********************************************************************//**
 * @brief GMWLObservation class extension
 ***************************************************************************/
%extend GMWLObservation {
    char *__str__() {
        static std::string result = self->print();
        return ((char*)result.c_str());
    }
    GMWLObservation copy() {
        return (*self);
    }
};
