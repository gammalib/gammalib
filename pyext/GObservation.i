/***************************************************************************
 *    GObservation.i  -  Observation abstract base class SWIG interface    *
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
 * @brief GObservation class SWIG file.
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
 * @brief Abstract interface for the observation classes.
 ***************************************************************************/
class GObservation {

    // Friend classes
    friend class GObservations;

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
    virtual double model(const GModels& models, const GEvent& event,
                         GVector* gradient = NULL) const;
    virtual double npred(const GModels& models, GVector* gradient = NULL) const;

    // Implemented methods
    void               obsname(const std::string& obsname);
    void               ebounds(const GEbounds& ebounds);
    void               gti(const GGti& gti);
    void               roi(const GRoi& roi);
    void               statistics(const std::string& statistics);
    const std::string& obsname(void) const;
    GTime              tstart(void) const;
    GTime              tstop(void) const;
    GEnergy            emin(void) const;
    GEnergy            emax(void) const;
    const GEbounds&    ebounds(void) const;
    const GGti&        gti(void) const;
    const GRoi*        roi(void) const;
    const GEvents*     events(void) const;
    const std::string& statistics(void) const;
};


/***********************************************************************//**
 * @brief GObservation class extension
 ***************************************************************************/
%extend GObservation {
    char *__str__() {
        return tochar(self->print());
    }
};
