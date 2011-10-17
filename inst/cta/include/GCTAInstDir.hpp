/***************************************************************************
 *            GCTAInstDir.hpp  -  CTA instrument direction class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *  This program is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU General Public License as published by   *
 *  the Free Software Foundation, either version 3 of the License, or      *
 *  (at your option) any later version.                                    *
 *                                                                         *
 *  This program is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *  GNU General Public License for more details.                           *
 *                                                                         *
 *  You should have received a copy of the GNU General Public License      *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GCTAInstDir.hpp
 * @brief CTA instrument direction class interface definition
 * @author J. Knodlseder
 */

#ifndef GCTAINSTDIR_HPP
#define GCTAINSTDIR_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GInstDir.hpp"
#include "GSkyDir.hpp"


/***********************************************************************//**
 * @class GCTAInstDir
 *
 * @brief CTA instrument direction class.
 *
 * The CTA instrument direction is an encapsulation of a sky direction
 * as CTA is an imaging device.
 ***************************************************************************/
class GCTAInstDir : public GInstDir {

public:
    // Constructors and destructors
    GCTAInstDir(void);
    explicit GCTAInstDir(const GSkyDir& dir);
    GCTAInstDir(const GCTAInstDir& dir);
    virtual ~GCTAInstDir(void);

    // Operators
    GCTAInstDir& operator= (const GCTAInstDir& dir);

    // Implemented pure virtual base class methods
    void         clear(void);
    GCTAInstDir* clone(void) const;
    std::string  print(void) const;

    // Other methods
    void         skydir(const GSkyDir& dir) { m_dir=dir; }
    void         radec(const double& ra, const double& dec) { m_dir.radec(ra,dec); }
    void         radec_deg(const double& ra, const double& dec) { m_dir.radec_deg(ra,dec); }
    void         lb(const double& l, const double& b) { m_dir.lb(l,b); }
    void         lb_deg(const double& l, const double& b) { m_dir.lb_deg(l,b); }
    void         rotate_deg(const double& phi, const double& theta);
    GSkyDir      skydir(void) const { return m_dir; }
    double       l(void) const { return m_dir.l(); }
    double       l_deg(void) const { return m_dir.l_deg(); }
    double       b(void) const { return m_dir.b(); }
    double       b_deg(void) const { return m_dir.b_deg(); }
    double       ra(void) const { return m_dir.ra(); }
    double       ra_deg(void) const { return m_dir.ra_deg(); }
    double       dec(void) const { return m_dir.dec(); }
    double       dec_deg(void) const { return m_dir.dec_deg(); }
    double       dist(const GSkyDir& dir) const { return m_dir.dist(dir); }
    double       dist_deg(const GSkyDir& dir) const { return m_dir.dist_deg(dir); }
    double       dist(const GCTAInstDir& dir) const;
    double       dist_deg(const GCTAInstDir& dir) const;
    double       posang(const GSkyDir& dir) const { return m_dir.posang(dir); }
    double       posang_deg(const GSkyDir& dir) const { return m_dir.posang_deg(dir); }
    double       posang(const GCTAInstDir& dir) const;
    double       posang_deg(const GCTAInstDir& dir) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCTAInstDir& dir);
    void free_members(void);

    // Data members
    GSkyDir m_dir;  //!< Observed incident direction of event
};

#endif /* GCTAINSTDIR_HPP */
