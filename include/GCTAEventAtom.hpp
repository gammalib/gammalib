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
#include "GEventAtom.hpp"
#include "GModels.hpp"
#include "GVector.hpp"


/***********************************************************************//**
 * @class GCTAEventAtom
 *
 * @brief GCTAEventAtom class interface defintion.
 ***************************************************************************/
class GCTAEventAtom : public GEventAtom {

    // Friend classes
    friend class GCTAEventList;

public:
    // Constructors and destructors
    GCTAEventAtom();
    GCTAEventAtom(const GCTAEventAtom& atom);
    virtual ~GCTAEventAtom();

    // Operators
    GCTAEventAtom& operator= (const GCTAEventAtom& atom);

    // Methods
    double model(GModels& models);
    double model(GModels& models, GVector* gradient);
    
protected:
    // Protected methods
    void           init_members(void);
    void           copy_members(const GCTAEventAtom& atom);
    void           free_members(void);
    GCTAEventAtom* clone(void) const;

    // CTA data format attributes
    long    m_event_id;             //!< Event identifier
    int     m_flags;                //!< Flags
    int     m_multip;               //!< Multiplicity
    int     m_telmask;              //!< 
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
    
    // Other attributes
};

#endif /* GCTAEVENTATOM_HPP */
