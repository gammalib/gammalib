/***************************************************************************
 *      GModelSpatialRadialShell.hpp - Radial shell source model class     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2022 by Christoph Deil                              *
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
 * @file GModelSpatialRadialShell.hpp
 * @brief Radial shell model class interface definition
 * @author Christoph Deil
 */

#ifndef GMODELSPATIALRADIALSHELL_HPP
#define GMODELSPATIALRADIALSHELL_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelSpatialRadial.hpp"
#include "GModelPar.hpp"

/* __ Forward declarations _______________________________________________ */
class GEnergy;
class GTime;
class GPhoton;
class GRan;
class GSkyDir;
class GSkyRegion;
class GXmlElement;


/**************************************************************************
 * @class GModelSpatialRadialShell
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
class GModelSpatialRadialShell : public GModelSpatialRadial {

public:
    // Constructors and destructors
    GModelSpatialRadialShell(void);
    GModelSpatialRadialShell(const bool& dummy, const std::string& type);
    GModelSpatialRadialShell(const GSkyDir&     dir,
                             const double&      radius,
                             const double&      width,
                             const std::string& coordsys = "CEL");
    explicit GModelSpatialRadialShell(const GXmlElement& xml);
    GModelSpatialRadialShell(const GModelSpatialRadialShell& model);
    virtual ~GModelSpatialRadialShell(void);

    // Operators
    virtual GModelSpatialRadialShell& operator=(const GModelSpatialRadialShell& model);

    // Implemented pure virtual methods
    virtual void                      clear(void);
    virtual GModelSpatialRadialShell* clone(void) const;
    virtual std::string               classname(void) const;
    virtual double                    eval(const double&  theta,
                                           const GEnergy& energy,
                                           const GTime&   time,
                                           const bool&    gradients = false) const;
    virtual GSkyDir                   mc(const GEnergy& energy,
                                         const GTime& time,
                                         GRan& ran) const;
    virtual bool                      contains(const GSkyDir& dir,
                                               const double&  margin = 0.0) const;
    virtual double                    theta_max(void) const;
    virtual void                      read(const GXmlElement& xml);
    virtual void                      write(GXmlElement& xml) const;
    virtual std::string               print(const GChatter& chatter = NORMAL) const;

    // Other methods
    double radius(void) const;
    double width(void) const;
    void   radius(const double& radius);
    void   width(const double& width);

protected:
    // Protected methods
    void          init_members(void);
    void          copy_members(const GModelSpatialRadialShell& model);
    void          free_members(void);
    void          update(void) const;
    static double f1(double x);
    static double f2(double x);
    virtual void  set_region(void) const;

    // Protected members
    GModelPar m_radius; //!< Inner shell radius (deg)
    GModelPar m_width;  //!< Shell thickness (deg)

    // Cached members used for pre-computations
    mutable double m_last_radius;   //!< Last shell radius (deg)
    mutable double m_last_width;    //!< Last shell width (deg)
    mutable double m_theta_in;      //!< Inner shell radius (rad)
    mutable double m_x_in;          //!< sin(m_theta_in)^2
    mutable double m_theta_out;     //!< Outer shell radius (rad)
    mutable double m_x_out;         //!< sin(m_theta_out)^2
    mutable double m_norm;          //!< Shell normalization
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GModelSpatialRadialShell").
 ***************************************************************************/
inline
std::string GModelSpatialRadialShell::classname(void) const
{
    return ("GModelSpatialRadialShell");
}


/***********************************************************************//**
 * @brief Return shell radius
 *
 * @return Shell radius (degrees).
 *
 * Returns the shell radius in degrees.
 ***************************************************************************/
inline
double GModelSpatialRadialShell::radius(void) const
{
    return (m_radius.value());
}


/***********************************************************************//**
 * @brief Set shell radius 
 *
 * @param[in] radius Shell radius (degrees).
 *
 * Sets the shell radius in degrees.
 ***************************************************************************/
inline
void GModelSpatialRadialShell::radius(const double& radius)
{
    m_radius.value(radius);
    return;
}


/***********************************************************************//**
 * @brief Return shell width
 *
 * @return Shell width (degrees).
 *
 * Returns the shell width in degrees.
 ***************************************************************************/
inline
double GModelSpatialRadialShell::width(void) const
{
    return (m_width.value());
}


/***********************************************************************//**
 * @brief Set width radius 
 *
 * @param[in] width Shell width (degrees).
 *
 * Sets the shell width in degrees.
 ***************************************************************************/
inline
void GModelSpatialRadialShell::width(const double& width)
{
    m_width.value(width);
    return;
}

#endif /* GMODELSPATIALRADIALSHELL_HPP */
