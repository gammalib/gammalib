/***************************************************************************
 *         GModelRadialShell.hpp - Radial shell source model class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2013 by Christoph Deil                              *
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
 * @file GModelRadialShell.hpp
 * @brief Radial shell model class interface definition
 * @author Christoph Deil
 */

#ifndef GMODELRADIALSHELL_HPP
#define GMODELRADIALSHELL_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelSpatialRadial.hpp"
#include "GModelPar.hpp"
#include "GSkyDir.hpp"
#include "GXmlElement.hpp"


/**************************************************************************
 * @class GModelRadialShell
 *
 * @brief Shell source model class
 *
 * This class implements the spatial component of the factorised source
 * model for a shell source (can be used e.g. as a toy supernova remnant
 * model). The shell is simply the volume between an inner and outer radius
 * (a large sphere with a smaller sphere cut out) with constant volume
 * emissivity and no absorption. To determine the surface brightness
 * distribution on the sky, the shell is analytically integrated along lines
 * of sight.
 * The shell is parametrised by the inner shell radius and the shell width.
 * Obviously, the shell width has to be a positive quantity. To assure
 * positivity, 1 arcsec is added internally to the shell width. This should
 * be negligible compared to shell widths encountered (and detectable) in
 * the gamma-ray domain, and improves the convergence of the fitting
 * algorithms.
 ***************************************************************************/
class GModelRadialShell : public GModelSpatialRadial {

public:
    // Constructors and destructors
    GModelRadialShell(void);
    explicit GModelRadialShell(const GSkyDir& dir,
                               const double& radius, const double& width,
                               const bool& small_angle = true);
    explicit GModelRadialShell(const GXmlElement& xml);
    GModelRadialShell(const GModelRadialShell& model);
    virtual ~GModelRadialShell(void);

    // Operators
    virtual GModelRadialShell& operator=(const GModelRadialShell& model);

    // Implemented pure virtual methods
    virtual void               clear(void);
    virtual GModelRadialShell* clone(void) const;
    virtual std::string        type(void) const { return "ShellFunction"; }
    virtual double             eval(const double& theta) const;
    virtual double             eval_gradients(const double& theta) const;
    virtual GSkyDir            mc(GRan& ran) const;
    virtual double             theta_max(void) const;
    virtual void               read(const GXmlElement& xml);
    virtual void               write(GXmlElement& xml) const;
    virtual std::string        print(void) const;

    // Other methods
    double  radius(void) const { return m_radius.real_value(); }
    double  width(void) const { return m_width.real_value(); }
    bool    small_angle(void) const { return m_small_angle; }
    void    radius(const double& radius) { m_radius.real_value(radius); }
    void    width(const double& width) { m_width.real_value(width); }
    void    small_angle(const bool& small_angle) { m_small_angle = small_angle; }

protected:
    // Protected methods
    void          init_members(void);
    void          copy_members(const GModelRadialShell& model);
    void          free_members(void);
    void          update(void) const;
    static double f1(double x);
    static double f2(double x);

    // Protected members
    GModelPar       m_radius;        //!< Inner shell radius (deg)
    GModelPar       m_width;         //!< Shell thickness (deg)
    bool            m_small_angle;   //!< Use small angle approximate formulae

    // Cached members used for pre-computations
    mutable double  m_last_radius;   //!< Last shell radius (deg)
    mutable double  m_last_width;    //!< Last shell width (deg)
    mutable double  m_theta_in;      //!< Inner shell radius (rad)
    mutable double  m_x_in;          //!< m_theta_in^2 or sin(m_theta_in)^2
    mutable double  m_theta_out;     //!< Outer shell radius (rad)
    mutable double  m_x_out;         //!< m_theta_out^2 or sin(m_theta_out)^2
    mutable double  m_norm;          //!< Shell normalization
};

#endif /* GMODELRADIALSHELL_HPP */
