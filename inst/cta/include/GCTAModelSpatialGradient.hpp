/***************************************************************************
 *      GCTAModelSpatialGradient.hpp - Spatial gradient CTA model class    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2018 by Juergen Knoedlseder                              *
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
 * @file GCTAModelSpatialGradient.hpp
 * @brief Spatial gradient CTA interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTAMODELSPATIALGRADIENT_HPP
#define GCTAMODELSPATIALGRADIENT_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelPar.hpp"
#include "GCTAModelSpatial.hpp"

/* __ Forward declarations _______________________________________________ */
class GRan;
class GXmlElement;
class GCTAInstDir;
class GEnergy;
class GTime;
class GCTAInstDir;
class GCTAObservation;


/***********************************************************************//**
 * @class GCTAModelSpatialGradient
 *
 * @brief Spatial gradient CTA model class
 ***************************************************************************/
class GCTAModelSpatialGradient  : public GCTAModelSpatial {

public:
    // Constructors and destructors
    GCTAModelSpatialGradient(void);
    GCTAModelSpatialGradient(const double& detx_gradient,
                             const double& dety_gradient);
    explicit GCTAModelSpatialGradient(const GXmlElement& xml);
    GCTAModelSpatialGradient(const GCTAModelSpatialGradient& model);
    virtual ~GCTAModelSpatialGradient(void);

    // Operators
    virtual GCTAModelSpatialGradient& operator=(const GCTAModelSpatialGradient& model);

    // Implemented pure virtual methods
    virtual void                      clear(void);
    virtual GCTAModelSpatialGradient* clone(void) const;
    virtual std::string               classname(void) const;
    virtual std::string               type(void) const;
    virtual double                    eval(const GCTAInstDir& dir,
                                           const GEnergy&     energy,
                                           const GTime&       time,
                                           const bool&        gradients = false) const;
    virtual double                    mc_max_value(const GCTAObservation& obs) const;
    virtual void                      read(const GXmlElement& xml);
    virtual void                      write(GXmlElement& xml) const;
    virtual std::string               print(const GChatter& chatter = NORMAL) const;

    // Other methods
    double detx_gradient(void) const;
    double dety_gradient(void) const;
    void   detx_gradient(const double& detx_gradient);
    void   dety_gradient(const double& dety_gradient);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCTAModelSpatialGradient& model);
    void free_members(void);

    // Protected members
    GModelPar m_detx_gradient;        //!< DETX gradient
    GModelPar m_dety_gradient;        //!< DETY gradient
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCTAModelSpatialGradient").
 ***************************************************************************/
inline
std::string GCTAModelSpatialGradient::classname(void) const
{
    return ("GCTAModelSpatialGradient");
}


/***********************************************************************//**
 * @brief Return model type
 *
 * @return Model type "Gradient".
 ***************************************************************************/
inline
std::string GCTAModelSpatialGradient::type(void) const
{
    return ("Gradient");
}


/***********************************************************************//**
 * @brief Return DETX gradient
 *
 * @return DETX gradient.
 ***************************************************************************/
inline
double GCTAModelSpatialGradient::detx_gradient(void) const
{
    return (m_detx_gradient.value());
}


/***********************************************************************//**
 * @brief Return DETY gradient
 *
 * @return DETY gradient.
 ***************************************************************************/
inline
double GCTAModelSpatialGradient::dety_gradient(void) const
{
    return (m_dety_gradient.value());
}


/***********************************************************************//**
 * @brief Set DETX gradient
 *
 * @param[in] detx_gradient DETX gradient.
 ***************************************************************************/
inline
void GCTAModelSpatialGradient::detx_gradient(const double& detx_gradient)
{
    m_detx_gradient.value(detx_gradient);
    return;
}


/***********************************************************************//**
 * @brief Set DETY gradient
 *
 * @param[in] dety_gradient DETY gradient.
 ***************************************************************************/
inline
void GCTAModelSpatialGradient::dety_gradient(const double& dety_gradient)
{
    m_dety_gradient.value(dety_gradient);
    return;
}

#endif /* GCTAMODELSPATIALGRADIENT_HPP */
