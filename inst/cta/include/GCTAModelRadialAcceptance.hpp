/***************************************************************************
 *      GCTAModelRadialAcceptance.hpp - Radial acceptance model class      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2018 by Juergen Knoedlseder                         *
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
 * @file GCTAModelRadialAcceptance.hpp
 * @brief Radial acceptance model class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTAMODELRADIALACCEPTANCE_HPP
#define GCTAMODELRADIALACCEPTANCE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <cmath>
#include "GModelData.hpp"
#include "GModelPar.hpp"
#include "GModelSpectral.hpp"
#include "GModelTemporal.hpp"
#include "GEvent.hpp"
#include "GObservation.hpp"
#include "GXmlElement.hpp"
#include "GFunction.hpp"
#include "GCTAEventList.hpp"
#include "GCTAModelRadial.hpp"


/***********************************************************************//**
 * @class GCTAModelRadialAcceptance
 *
 * @brief Radial acceptance model class
 *
 * This class implements a radial acceptance model for CTA.
 ***************************************************************************/
class GCTAModelRadialAcceptance : public GModelData {

public:
    // Constructors and destructors
    GCTAModelRadialAcceptance(void);
    explicit GCTAModelRadialAcceptance(const GXmlElement& xml);
    GCTAModelRadialAcceptance(const GCTAModelRadial& radial,
                              const GModelSpectral&  spectral);
    GCTAModelRadialAcceptance(const GCTAModelRadial& radial,
                              const GModelSpectral&  spectral,
                              const GModelTemporal&  temporal);
    GCTAModelRadialAcceptance(const GCTAModelRadialAcceptance& model);
    virtual ~GCTAModelRadialAcceptance(void);

    // Operators
    virtual GCTAModelRadialAcceptance& operator=(const GCTAModelRadialAcceptance& model);

    // Implemented pure virtual methods
    virtual void                       clear(void);
    virtual GCTAModelRadialAcceptance* clone(void) const;
    virtual std::string                classname(void) const;
    virtual std::string                type(void) const;
    virtual bool                       is_constant(void) const;
    virtual double                     eval(const GEvent& event,
                                            const GObservation& obs,
                                            const bool& gradients = false) const;
    virtual double                     npred(const GEnergy& obsEng, const GTime& obsTime,
                                             const GObservation& obs) const;
    virtual GCTAEventList*             mc(const GObservation& obs, GRan& ran) const;
    virtual void                       read(const GXmlElement& xml);
    virtual void                       write(GXmlElement& xml) const;
    virtual std::string                print(const GChatter& chatter = NORMAL) const;

    // Other methods
    GCTAModelRadial* radial(void)   const;
    GModelSpectral*  spectral(void) const;
    GModelTemporal*  temporal(void) const;
    void             radial(const GCTAModelRadial* radial);
    void             spectral(const GModelSpectral* spectral);
    void             temporal(const GModelTemporal* temporal);

protected:
    // Protected methods
    void             init_members(void);
    void             copy_members(const GCTAModelRadialAcceptance& model);
    void             free_members(void);
    void             set_pointers(void);
    bool             valid_model(void) const;
    GCTAModelRadial* xml_radial(const GXmlElement& radial) const;
    GModelSpectral*  xml_spectral(const GXmlElement& spectral) const;
    GModelTemporal*  xml_temporal(const GXmlElement& temporal) const;

    // ROI integration kernel
    class roi_kern : public GFunction {
    public:
        roi_kern(const GCTAModelRadial* parent, const double& roi, const double& dist) :
                 m_parent(parent),
                 m_roi(roi),
                 m_cosroi(std::cos(roi)),
                 m_dist(dist),
                 m_cosdist(std::cos(dist)),
                 m_sindist(std::sin(dist)) { }
        double eval(const double& r);
    protected:
        const GCTAModelRadial* m_parent;   //!< Pointer to radial model
        double                 m_roi;      //!< ROI radius in radians
        double                 m_cosroi;   //!< Cosine of ROI radius
        double                 m_dist;     //!< Distance between pointing and ROI centre in radians
        double                 m_cosdist;  //!< Cosine of distance
        double                 m_sindist;  //!< Sinus of distance
    };

    // Proteced data members
    GCTAModelRadial* m_radial;       //!< Radial model
    GModelSpectral*  m_spectral;     //!< Spectral model
    GModelTemporal*  m_temporal;     //!< Temporal model
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCTAModelRadialAcceptance").
 ***************************************************************************/
inline
std::string GCTAModelRadialAcceptance::classname(void) const
{
    return ("GCTAModelRadialAcceptance");
}


/***********************************************************************//**
 * @brief Return model type
 *
 * @return Model type.
 *
 * Returns the type of the model. The type for a radial acceptance model is
 * "RadialAcceptance".
 ***************************************************************************/
inline
std::string GCTAModelRadialAcceptance::type(void) const
{
    return ("RadialAcceptance");
}


/***********************************************************************//**
 * @brief Signals if model is temporally constant
 *
 * @return True if model is temporally constant, false otherwise.
 *
 * Signals if the model is temporally constant. A temporally constant model
 * is a model that has a temporal component of type "Constant".
 ***************************************************************************/
inline
bool GCTAModelRadialAcceptance::is_constant(void) const
{
    return (m_temporal != NULL && m_temporal->type() == "Constant");
}


/***********************************************************************//**
 * @brief Return radial model component
 *
 * @return Pointer to radial model component.
 *
 * Returns a pointer to the radial model component of the model. The pointer
 * is of type GCTAModelRadial. Note that a NULL pointer may be returned if
 * the model has no radial model component.
 ***************************************************************************/
inline
GCTAModelRadial* GCTAModelRadialAcceptance::radial(void) const
{
    return (m_radial);
}


/***********************************************************************//**
 * @brief Return spectral model component
 *
 * @return Pointer to spectral model component.
 *
 * Returns a pointer to the spectral model component of the model. The
 * pointer is of type GModelSpectral. Note that a NULL pointer may be
 * returned if the model has no spectral model component.
 ***************************************************************************/
inline
GModelSpectral* GCTAModelRadialAcceptance::spectral(void) const
{
    return (m_spectral);
}


/***********************************************************************//**
 * @brief Return temporal model component
 *
 * @return Pointer to temporal model component.
 *
 * Returns a pointer to the temporal model component of the model. The
 * pointer is of type GModelTemporal. Note that a NULL pointer may be
 * returned if the model has no temporal model component.
 ***************************************************************************/
inline
GModelTemporal* GCTAModelRadialAcceptance::temporal(void) const
{
    return (m_temporal);
}

#endif /* GCTAMODELRADIALACCEPTANCE_HPP */
