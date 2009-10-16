/***************************************************************************
 *                GLATEventAtom.hpp  -  LAT event atom class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GLATEventAtom.hpp
 * @brief GLATEventAtom class interface definition.
 * @author J. Knodlseder
 */

#ifndef GLATEVENTATOM_HPP
#define GLATEVENTATOM_HPP

/* __ Includes ___________________________________________________________ */
#include "GEventAtom.hpp"
#include "GModels.hpp"
#include "GVector.hpp"


/***********************************************************************//**
 * @class GLATEventAtom
 *
 * @brief GLATEventAtom class interface defintion.
 ***************************************************************************/
class GLATEventAtom : public GEventAtom {

    // Friend classes
    friend class GLATEventList;

public:
    // Constructors and destructors
    GLATEventAtom();
    GLATEventAtom(const GLATEventAtom& atom);
    virtual ~GLATEventAtom();

    // Operators
    GLATEventAtom& operator= (const GLATEventAtom& atom);

    // Methods
    double model(GModels& models, GVector* gradient) const;
    
protected:
    // Protected methods
    void           init_members(void);
    void           copy_members(const GLATEventAtom& atom);
    void           free_members(void);
    GLATEventAtom* clone(void) const;

    // LAT data format attributes
    float   m_theta;                //!< Zenith angle in instrument system
    float   m_phi;                  //!< Azimuth angle in instrument system
    float   m_zenith_angle;         //!< Zenith angle in Earth system
    float   m_earth_azimuth_angle;  //!< Azimuth angle in Earth system
    long    m_event_id;             //!< ID number of original event
    long    m_run_id;               //!< Run number of original event
    short   m_recon_version;        //!< Version of event reconstruction software
    short   m_calib_version[3];     //!< Version of calibration tables for ACD, CAL
    short   m_event_class;          //!< Event class: 0, 1, 2, ...
    short   m_conversion_type;      //!< Type of conversion: 0=Front, 1=Back
    double  m_livetime;             //!< Accumulated livetime since mission start
    double* m_difrsp;               //!< Diffuse response components
    
    // Other attributes
    int     m_num_difrsp;           //!< Number of diffuse model components
};

#endif /* GLATEVENTATOM_HPP */
