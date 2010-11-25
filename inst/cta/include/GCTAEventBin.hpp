/***************************************************************************
 *                 GCTAEventBin.hpp  -  CTA event bin class                *
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
 * @file GCTAEventBin.hpp
 * @brief GCTAEventBin class interface definition.
 * @author J. Knodlseder
 */

#ifndef GCTAEVENTBIN_HPP
#define GCTAEVENTBIN_HPP

/* __ Includes ___________________________________________________________ */
#include <ostream>
#include "GEventBin.hpp"
#include "GModels.hpp"
#include "GVector.hpp"
#include "GEnergy.hpp"
#include "GCTAInstDir.hpp"
#include "GCTAPointing.hpp"
#include "GCTAResponse.hpp"


/***********************************************************************//**
 * @class GCTAEventBin
 *
 * @brief GCTAEventBin class interface defintion.
 *
 * The CTA event bin class holds all information that is relevant for the
 * analysis of an event cube bin. It derives from the abstract GEventBin
 * base class that holds pointers to the actual number of counts, to the
 * mean time and to the mean energy of the specific bin. These pointers
 * are set by the GCTAEventCube class that is the container class for
 * CTA events. Note that the counts are in fact stored in GCTAEventCube.
 * GCTAEventBin just provides an interface to access event information bin
 * by bin.
 ***************************************************************************/
class GCTAEventBin : public GEventBin {

    // Friend classes
    friend class GCTAEventCube;

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GCTAEventBin& bin);

public:
    // Constructors and destructors
    GCTAEventBin(void);
    GCTAEventBin(const GCTAEventBin& bin);
    virtual ~GCTAEventBin(void);

    // Operators
    GCTAEventBin& operator= (const GCTAEventBin& bin);

    // Methods
    double              model(GModels& models, GVector* gradient = NULL) const;
    double              size(void) const;
    const GCTAInstDir*  dir(void) const { return m_dir; }
    const GCTAPointing* pnt(void) const { return m_pnt; }
    const GCTAResponse* rsp(void) const { return m_rsp; }
    const double*       omega(void) const { return m_omega; }
    const GEnergy*      ewidth(void) const { return m_ewidth; }
    const double*       ontime(void) const { return m_ontime; }
    
protected:
    // Protected methods
    void          init_members(void);
    void          copy_members(const GCTAEventBin& bin);
    void          free_members(void);
    GCTAEventBin* clone(void) const;

    // CTA specific event attributes
    GCTAInstDir*  m_dir;     //!< Pointer to event direction
    GCTAPointing* m_pnt;     //!< Pointer to instrument pointing
    GCTAResponse* m_rsp;     //!< Pointer to instrument response function
    double*       m_omega;   //!< Pointer to solid angle of pixel (sr)
    GEnergy*      m_ewidth;  //!< Pointer to energy width of bin
    double*       m_ontime;  //!< Pointer to ontime of bin (seconds)

};

#endif /* GCTAEVENTBIN_HPP */
