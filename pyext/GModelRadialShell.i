/***************************************************************************
 *     GModelRadialShell.i  -  Radial spatial shell source model class     *
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
 * @file GModelRadialShell.i
 * @brief Radial spatial shell model class Python interface definition
 * @author C. Deil
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModelRadialShell.hpp"
#include "GTools.hpp"
%}

/**************************************************************************
 * @class GModelRadialShell
 *
 * @brief Radial shell source model class
 ***************************************************************************/
class GModelRadialShell : public GModelRadial {
public:
    // Constructors and destructors
    GModelRadialShell(void);
    explicit GModelRadialShell(const GSkyDir& dir,
                               const double& radius, const double& width,
                               const bool& small_angle = true);
    explicit GModelRadialShell(const GXmlElement& xml);
    GModelRadialShell(const GModelRadialShell& model);
    virtual ~GModelRadialShell(void);

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

    // Other methods
    double  radius(void) const { return m_radius.real_value(); }
    double  width(void) const { return m_width.real_value(); }
    bool    small_angle(void) const { return m_small_angle; }
    void    radius(const double& radius) { m_radius.real_value(radius); }
    void    width(const double& width) { m_width.real_value(width); }
    void    small_angle(const bool& small_angle) { m_small_angle = small_angle; }
};


/***********************************************************************//**
 * @brief GModelRadialShell class extension
 ***************************************************************************/
%extend GModelRadialShell {
    GModelRadialShell copy() {
        return (*self);
    }
};


/***********************************************************************//**
 * @brief GModelRadialShell type casts
 ***************************************************************************/
%inline %{
    GModelRadialShell* cast_GModelRadialShell(GModelSpatial* model) {
        return dynamic_cast<GModelRadialShell*>(model);
    }
%};
