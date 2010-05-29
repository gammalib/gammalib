/***************************************************************************
 *                GCTAEventAtom.hpp  -  CTA event atom class               *
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
 * @file GCTAEventAtom.hpp
 * @brief GCTAEventAtom class interface definition.
 * @author J. Knodlseder
 */

#ifndef GCTAEVENTATOM_HPP
#define GCTAEVENTATOM_HPP

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include "GEventAtom.hpp"
#include "GModels.hpp"
#include "GVector.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GCTAInstDir.hpp"
#include "GCTAPointing.hpp"
#include "GCTAResponse.hpp"


/***********************************************************************//**
 * @class GCTAEventAtom
 *
 * @brief GCTAEventAtom class interface defintion.
 ***************************************************************************/
class GCTAEventAtom : public GEventAtom {

    // Friend classes
    friend class GCTAEventList;

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GCTAEventAtom& atom);

public:
    // Constructors and destructors
    GCTAEventAtom(void);
    GCTAEventAtom(const GCTAEventAtom& atom);
    virtual ~GCTAEventAtom(void);

    // Operators
    GCTAEventAtom& operator= (const GCTAEventAtom& atom);

    // Methods
    double              model(GModels& models) const;
    double              model(GModels& models, GVector* gradient) const;
    const GCTAInstDir*  dir(void) const { return &m_dir; }
    const GCTAPointing* pnt(void) const { return &m_pnt; }
    const GCTAResponse* rsp(void) const { return m_rsp; }
    
protected:
    // Protected methods
    void           init_members(void);
    void           copy_members(const GCTAEventAtom& atom);
    void           free_members(void);
    GCTAEventAtom* clone(void) const;
    
    // CTA attributes
    GCTAInstDir   m_dir;            //!< Event direction
    GCTAPointing  m_pnt;            //!< CTA pointing
    GCTAResponse* m_rsp;            //!< Pointer to CTA instrument response function

    // CTA data format attributes
    long    m_event_id;             //!< Event identifier
    char    m_flags;                //!< Flags
    int     m_multip;               //!< Multiplicity
    char    m_telmask;              //!< Telescope mask
    float   m_dir_err;              //!< Error on event direction
    float   m_detx;                 //!<
    float   m_dety;                 //!<
    float   m_alt_pnt;              //!< Pointing altitude
    float   m_az_pnt;               //!< Pointing azimuth
    float   m_alt;                  //!< Event altitude
    float   m_az;                   //!< Event azimuth
    float   m_corex;                //!< 
    float   m_corey;                //!<
    float   m_core_err;             //!<
    float   m_xmax;                 //!<
    float   m_xmax_err;             //!<
    float   m_energy_err;           //!< Error on event energy (MeV)
    float   m_hil_msw;              //!< Hillas width
    float   m_hil_msw_err;          //!< Error on Hillas width
    float   m_hil_msl;              //!< Hillas length
    float   m_hil_msl_err;          //!< Error on Hillas length
};

#endif /* GCTAEVENTATOM_HPP */
