/***************************************************************************
 *               GResponse.i  -  Abstract response base class              *
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
 * @file GResponse.i
 * @brief Abstract response base class python bindings
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GResponse.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GResponse
 *
 * @brief Abstract response base class python bindings
 ***************************************************************************/
class GResponse {
public:
    // Constructors and destructors
    GResponse(void);
    GResponse(const GResponse& rsp);
    virtual ~GResponse(void);

    // Pure virtual methods
    virtual void       clear(void) = 0;
    virtual GResponse* clone(void) const = 0;
    virtual bool       hasedisp(void) const = 0;
    virtual bool       hastdisp(void) const = 0;
    virtual double     irf(const GEvent& event, const GSkyModel& model,
                           const GEnergy& srcEng, const GTime& srcTime,
                           const GObservation& obs) const = 0;
    virtual double     npred(const GSkyModel& model, const GEnergy& srcEng,
                             const GTime& srcTime,
                             const GObservation& obs) const = 0;
};


/***********************************************************************//**
 * @brief GResponse class extension
 ***************************************************************************/
%extend GResponse {
    char *__str__() {
        return tochar(self->print());
    }
};
