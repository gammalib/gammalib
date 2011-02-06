/***************************************************************************
 *               GCTAObservation.hpp  -  CTA Observation class             *
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
 * @file GCTAObservation.hpp
 * @brief GCTAObservation class interface definition.
 * @author J. Knodlseder
 */

#ifndef GCTAOBSERVATION_HPP
#define GCTAOBSERVATION_HPP

/* __ Includes ___________________________________________________________ */
#include "GObservation.hpp"
#include "GCTAPointing.hpp"
#include "GCTAResponse.hpp"
#include "GTime.hpp"
#include "GModel.hpp"
#include "GFitsTable.hpp"


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

    // Operators
    GCTAObservation& operator= (const GCTAObservation& obs);

    // Implemented pure virtual base class methods
    void             clear(void);
    GCTAObservation* clone(void) const;
    GCTAResponse*    response(void) const;
    GCTAPointing*    pointing(const GTime& time) const;
    std::string      instrument(void) const;
    std::string      print(void) const;

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

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCTAObservation& obs);
    void free_members(void);
    void read_attributes(GFitsHDU* hdu);
    void read_ds_ebounds(GFitsHDU* hdu);
    void read_ds_roi(GFitsHDU* hdu);
    void write_attributes(GFitsHDU* hdu) const;

    // Npred integration methods
    double npred_temp(const GModel& model) const;

    // Protected members
    GCTAResponse* m_response;   //!< Pointer to instrument response functions
    GCTAPointing* m_pointing;   //!< Pointer to pointing direction
    int           m_obs_id;     //!< Observation ID
    double        m_livetime;   //!< Livetime
    double        m_ra_obj;     //!< Right Ascension of object
    double        m_dec_obj;    //!< Declination of object
};

#endif /* GCTAOBSERVATION_HPP */
