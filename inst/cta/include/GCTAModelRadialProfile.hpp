/***************************************************************************
 *       GCTAModelRadialProfile.hpp - Radial Profile CTA model class       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2013 by Juergen Knoedlseder                         *
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
 * @file GCTAModelRadialProfile.hpp
 * @brief Radial Profile model class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTAMODELRADIALPROFILE_HPP
#define GCTAMODELRADIALPROFILE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <cmath>
#include "GMath.hpp"
#include "GModelPar.hpp"
#include "GXmlElement.hpp"
#include "GCTAModelRadial.hpp"
#include "GCTAInstDir.hpp"
#include "GFunction.hpp"


/***********************************************************************//**
 * @class GCTAModelRadialProfile
 *
 * @brief Radial Profile CTA model class
 *
 * This class implements the radial profile function
 * \f[f(\theta) = (1 + (\theta/c_0)^{c_1})^{-c_2/c_1}\f]
 * where
 * \f$\theta\f$ is the offset angle (in degrees),
 * \f$c_0\f$ is the width of the profile (width),
 * \f$c_1\f$ is the width of the central plateau (core), and
 * \f$c_2\f$ is the size of the tail (tail).
 * Note that all 3 parameters are positive values.
 ***************************************************************************/
class GCTAModelRadialProfile  : public GCTAModelRadial {

public:
    // Constructors and destructors
    GCTAModelRadialProfile(void);
    explicit GCTAModelRadialProfile(const double& width, const double& core,
                                    const double& tail);
    explicit GCTAModelRadialProfile(const GXmlElement& xml);
    GCTAModelRadialProfile(const GCTAModelRadialProfile& model);
    virtual ~GCTAModelRadialProfile(void);

    // Operators
    virtual GCTAModelRadialProfile& operator= (const GCTAModelRadialProfile& model);

    // Implemented pure virtual methods
    virtual void                    clear(void);
    virtual GCTAModelRadialProfile* clone(void) const;
    virtual std::string             type(void) const;
    virtual double                  eval(const double& offset) const;
    virtual double                  eval_gradients(const double& offset) const;
    virtual GCTAInstDir             mc(const GCTAInstDir& dir, GRan& ran) const;
    virtual double                  omega(void) const;
    virtual void                    read(const GXmlElement& xml);
    virtual void                    write(GXmlElement& xml) const;
    virtual std::string             print(const GChatter& chatter = NORMAL) const;

    // Other methods
    double width(void) const;
    double core(void) const;
    double tail(void) const;
    void   width(const double& width);
    void   core(const double& core);
    void   tail(const double& tail);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCTAModelRadialProfile& model);
    void free_members(void);

    // Radial integration class (used by omega() method). Note that the
    // integration is done in radians
    class integrand : public GFunction {
    public:
        integrand(const GCTAModelRadialProfile* model) : m_model(model) { }
        double eval(const double& x) {
            return (std::sin(x)*m_model->eval(x*gammalib::rad2deg));
        }
    private:
        const GCTAModelRadialProfile* m_model;
    };

    // Protected members
    GModelPar m_width;        //!< Width parameter
    GModelPar m_core;         //!< Core parameter
    GModelPar m_tail;         //!< Tail parameter
};


/***********************************************************************//**
 * @brief Return model type
 *
 * @return Model type "Profile".
 ***************************************************************************/
inline
std::string GCTAModelRadialProfile::type(void) const
{
    return ("Profile");
}


/***********************************************************************//**
 * @brief Return profile width
 *
 * @return Profile width.
 ***************************************************************************/
inline
double GCTAModelRadialProfile::width(void) const
{
    return (m_width.value());
}


/***********************************************************************//**
 * @brief Return profile core
 *
 * @return Profile core.
 ***************************************************************************/
inline
double GCTAModelRadialProfile::core(void) const
{
    return (m_core.value());
}


/***********************************************************************//**
 * @brief Return profile tail
 *
 * @return Profile tail.
 ***************************************************************************/
inline
double GCTAModelRadialProfile::tail(void) const
{
    return (m_tail.value());
}


/***********************************************************************//**
 * @brief Set profile width
 *
 * @param[in] Profile width.
 ***************************************************************************/
inline
void GCTAModelRadialProfile::width(const double& width)
{
    m_width.value(width);
    return;
}


/***********************************************************************//**
 * @brief Set profile core
 *
 * @param[in] Profile core.
 ***************************************************************************/
inline
void GCTAModelRadialProfile::core(const double& core)
{
    m_core.value(core);
    return;
}


/***********************************************************************//**
 * @brief Set profile tail
 *
 * @param[in] Profile tail.
 ***************************************************************************/
inline
void GCTAModelRadialProfile::tail(const double& tail)
{
    m_tail.value(tail);
    return;
}

#endif /* GCTAMODELRADIALPROFILE_HPP */
