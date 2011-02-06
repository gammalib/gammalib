/***************************************************************************
 *       GCTAObservation.i  -  CTA Observation class SWIG interface        *
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
 * @file GCTAObservation.i
 * @brief GCTAObservation class SWIG interface.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAObservation.hpp"
%}


/***********************************************************************//**
 * @class GCTAObservation
 *
 * @brief Interface class for CTA observations.
 ***************************************************************************/
class GCTAObservation : public GObservation {
public:
    // Constructors and destructors
    GCTAObservation(void);
    GCTAObservation(const GCTAObservation& obs);
    virtual ~GCTAObservation(void);

    // Pure virtual base class methods
    void             clear(void);
    GCTAObservation* clone(void) const;
    GCTAResponse*    response(void) const;
    GCTAPointing*    pointing(const GTime& time) const;
    std::string      instrument(void) const;

    // Other methods
    void   load_unbinned(const std::string& filename);
    void   load_binned(const std::string& filename);
    void   save(const std::string& filename, bool clobber) const;
    void   response(const std::string& irfname, std::string caldb = "");
    void   pointing(const GCTAPointing& pointing);
    void   obs_id(const int& id) { m_obs_id=id; }
    void   ra_obj(const double& ra) { m_ra_obj=ra; }
    void   dec_obj(const double& dec) { m_dec_obj=dec; }
    int    obs_id(void) const { return m_obs_id; }
    double livetime(void) const { return m_livetime; }
    double ra_obj(void) const { return m_ra_obj; }
    double dec_obj(void) const { return m_dec_obj; }
};


/***********************************************************************//**
 * @brief GCTAObservation class extension
 ***************************************************************************/
%extend GCTAObservation {
    GCTAObservation copy() {
        return (*self);
    }
};
