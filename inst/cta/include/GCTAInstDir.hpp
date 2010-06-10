/***************************************************************************
 *            GCTAInstDir.hpp  -  CTA instrument direction class           *
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
 * @file GCTAInstDir.hpp
 * @brief GCTAInstDir class definition.
 * @author J. Knodlseder
 */

#ifndef GCTAINSTDIR_HPP
#define GCTAINSTDIR_HPP

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include "GInstDir.hpp"
#include "GSkyDir.hpp"


/***********************************************************************//**
 * @class GCTAInstDir
 *
 * @brief Interface for the CTA instrument direction class.
 *
 * The CTA instrument direction is an encapsulation of a sky direction
 * as CTA is an imaging device.
 ***************************************************************************/
class GCTAInstDir : public GInstDir {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GCTAInstDir& dir);

public:
    // Constructors and destructors
    GCTAInstDir(void);
    explicit GCTAInstDir(const GSkyDir& dir);
    GCTAInstDir(const GCTAInstDir& dir);
    virtual ~GCTAInstDir(void);

    // Operators
    GCTAInstDir& operator= (const GCTAInstDir& dir);

    // Methods
    void   clear(void);
    void   radec(const double& ra, const double& dec) { m_dir.radec(ra,dec); }
    void   radec_deg(const double& ra, const double& dec) { m_dir.radec_deg(ra,dec); }
    void   lb(const double& l, const double& b) { m_dir.lb(l,b); }
    void   lb_deg(const double& l, const double& b) { m_dir.lb_deg(l,b); }
    double l(void) const { return m_dir.l(); }
    double l_deg(void) const { return m_dir.l_deg(); }
    double b(void) const { return m_dir.b(); }
    double b_deg(void) const { return m_dir.b_deg(); }
    double ra(void) const { return m_dir.ra(); }
    double ra_deg(void) const { return m_dir.ra_deg(); }
    double dec(void) const { return m_dir.dec(); }
    double dec_deg(void) const { return m_dir.dec_deg(); }
    double dist(const GSkyDir& dir) const { return m_dir.dist(dir); }
    double dist_deg(const GSkyDir& dir) const { return m_dir.dist_deg(dir); }
    double dist(const GCTAInstDir& dir) const;
    double dist_deg(const GCTAInstDir& dir) const;

protected:
    // Protected methods
    void                 init_members(void);
    void                 copy_members(const GCTAInstDir& dir);
    void                 free_members(void);
    virtual GCTAInstDir* clone(void) const;
    
    // Data members
    GSkyDir m_dir;  //!< Incident direction of event
};

#endif /* GCTAINSTDIR_HPP */
