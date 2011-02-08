/***************************************************************************
 *       GObservation.i  -  Abstract observation abstract base class       *
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
 * @file GObservation.i
 * @brief Abstract observation base class Python interface definition
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GObservation.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GObservation
 *
 * @brief Abstract Python interface for the observation classes
 ***************************************************************************/
class GObservation {
public:
    // Constructors and destructors
    GObservation(void);
    GObservation(const GObservation& obs);
    virtual ~GObservation(void);

    // Pure virtual methods
    virtual void          clear(void) = 0;
    virtual GObservation* clone(void) const = 0;
    virtual GResponse*    response(void) const = 0;
    virtual GPointing*    pointing(const GTime& time) const = 0;
    virtual std::string   instrument(void) const = 0;

    // Virtual methods
    virtual double        model(const GModels& models, const GEvent& event,
                                GVector* gradient = NULL) const;
    virtual double        npred(const GModels& models, GVector* gradient = NULL) const;

    // Implemented methods
    void                  name(const std::string& name);
    void                  events(const GEvents* events);
    void                  statistics(const std::string& statistics);
    const std::string&    name(void) const;
    const GEvents*        events(void) const;
    const std::string&    statistics(void) const;
};


/***********************************************************************//**
 * @brief GObservation class extension
 ***************************************************************************/
%extend GObservation {
    char *__str__() {
        return tochar(self->print());
    }
};
