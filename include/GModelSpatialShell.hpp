/***************************************************************************
 *     GModelSpatialShell.hpp  -  Spatial shell source model class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011 by Christoph Deil                                   *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GModelSpatialShell.hpp
 * @brief Spatial shell model class interface definition
 * @author C. Deil
 */

#ifndef GMODELSPATIALSHELL_HPP
#define GMODELSPATIALSHELL_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelSpatial.hpp"
#include "GModelPar.hpp"
#include "GSkyDir.hpp"
#include "GXmlElement.hpp"


/**************************************************************************
 * @class GModelSpatialShell
 *
 * @brief Shell source model class
 *
 * This class implements the spatial component of the factorised source
 * model for a shell source (can be used e.g. as a toy supernova remnant model).
 * The shell is simply the volume between an inner and outer radius
 * (a large sphere with a smaller sphere cut out) with constant
 * volume emissivity and no absorption.
 * To get the surface brightness distribution on the sky, this sphere is
 * integrated along parallel lines of sight.
 ***************************************************************************/
class GModelSpatialShell : public GModelSpatial {

public:
    // Constructors and destructors
    GModelSpatialShell(void);
    explicit GModelSpatialShell(const GSkyDir& dir,
                                const double& radius, const double& width,
                                const bool& small_angle = true);
    explicit GModelSpatialShell(const GXmlElement& xml);
    GModelSpatialShell(const GModelSpatialShell& model);
    virtual ~GModelSpatialShell(void);

    // Operators
    virtual GModelSpatialShell& operator=(const GModelSpatialShell& model);

    // Implemented pure virtual methods
    virtual void                clear(void);
    virtual GModelSpatialShell* clone(void) const;
    virtual std::string         type(void) const { return "ShellFunction"; }
    virtual double              eval(const GSkyDir& srcDir) const;
    virtual double              eval_gradients(const GSkyDir& srcDir) const;
    virtual GSkyDir             mc(GRan& ran) const;
    virtual void                read(const GXmlElement& xml);
    virtual void                write(GXmlElement& xml) const;
    virtual std::string         print(void) const;

    // Other methods
    double  ra(void) const { return m_ra.real_value(); }
    double  dec(void) const { return m_dec.real_value(); }
    double  radius(void) const { return m_radius.real_value(); }
    double  width(void) const { return m_width.real_value(); }
    bool    small_angle(void) const { return m_small_angle; }
    GSkyDir dir(void) const;
    void    dir(const GSkyDir& dir);
    void    radius(const double& radius) { m_radius.real_value(radius); }
    void    width(const double& width) { m_width.real_value(width); }
    void    small_angle(const bool& small_angle) { m_small_angle = small_angle; }

protected:
    // Protected methods
    void          init_members(void);
    void          copy_members(const GModelSpatialShell& model);
    void          free_members(void);
    void          update() const;
    static double f1(double x);
    static double f2(double x);    

    // Protected members
    GModelPar       m_ra;            //!< Right Ascension of shell centre (deg)
    GModelPar       m_dec;           //!< Declination of shell centre (deg)
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

#endif /* GMODELSPATIALSHELL_HPP */
