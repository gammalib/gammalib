/***************************************************************************
 *           GCTAModelBackground.hpp - Background model class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2018-2021 by Juergen Knoedlseder                         *
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
 * @file GCTAModelBackground.hpp
 * @brief Background model class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTAMODELBACKGROUND_HPP
#define GCTAMODELBACKGROUND_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <cmath>
#include "GModelData.hpp"
#include "GModelSpectral.hpp"
#include "GModelTemporal.hpp"
#include "GFunction.hpp"
#include "GCTAModelSpatial.hpp"
#include "GCTAEventList.hpp"

/* __ Forward declarations _______________________________________________ */
class GEvent;
class GObservation;
class GModelPar;
class GXmlElement;


/***********************************************************************//**
 * @class GCTAModelBackground
 *
 * @brief Background model class
 *
 * This class implements a background model for CTA.
 ***************************************************************************/
class GCTAModelBackground : public GModelData {

public:
    // Constructors and destructors
    GCTAModelBackground(void);
    explicit GCTAModelBackground(const GXmlElement& xml);
    GCTAModelBackground(const GCTAModelSpatial& spatial,
                        const GModelSpectral&   spectral);
    GCTAModelBackground(const GCTAModelSpatial& spatial,
                        const GModelSpectral&   spectral,
                        const GModelTemporal&   temporal);
    GCTAModelBackground(const GCTAModelBackground& model);
    virtual ~GCTAModelBackground(void);

    // Operators
    virtual GCTAModelBackground& operator=(const GCTAModelBackground& model);

    // Implemented pure virtual methods
    virtual void                 clear(void);
    virtual GCTAModelBackground* clone(void) const;
    virtual std::string          classname(void) const;
    virtual std::string          type(void) const;
    virtual bool                 is_constant(void) const;
    virtual double               eval(const GEvent&       event,
                                      const GObservation& obs,
                                      const bool&         gradients = false) const;
    virtual double               npred(const GEnergy&       energy,
                                       const GTime&         time,
                                       const GPolarization& polarization,
                                       const GObservation&  obs) const;
    virtual GCTAEventList*       mc(const GObservation& obs, GRan& ran) const;
    virtual void                 read(const GXmlElement& xml);
    virtual void                 write(GXmlElement& xml) const;
    virtual std::string          print(const GChatter& chatter = NORMAL) const;

    // Other methods
    GCTAModelSpatial* spatial(void) const;
    GModelSpectral*   spectral(void) const;
    GModelTemporal*   temporal(void) const;
    void              spatial(const GCTAModelSpatial* spatial);
    void              spectral(const GModelSpectral* spectral);
    void              temporal(const GModelTemporal* temporal);

protected:
    // Protected methods
    void              init_members(void);
    void              copy_members(const GCTAModelBackground& model);
    void              free_members(void);
    void              set_pointers(void);
    GCTAModelSpatial* xml_spatial(const GXmlElement& spatial) const;
    GModelSpectral*   xml_spectral(const GXmlElement& spectral) const;
    GModelTemporal*   xml_temporal(const GXmlElement& temporal) const;
    bool              valid_model(void) const;

    // Proteced data members
    GCTAModelSpatial* m_spatial;      //!< Spatial model
    GModelSpectral*   m_spectral;     //!< Spectral model
    GModelTemporal*   m_temporal;     //!< Temporal model
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCTAModelBackground").
 ***************************************************************************/
inline
std::string GCTAModelBackground::classname(void) const
{
    return ("GCTAModelBackground");
}


/***********************************************************************//**
 * @brief Return model type
 *
 * @return Model type.
 *
 * Returns "CTABackground" as the type of the model.
 ***************************************************************************/
inline
std::string GCTAModelBackground::type(void) const
{
    return ("CTABackground");
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
bool GCTAModelBackground::is_constant(void) const
{
    return (m_temporal != NULL && m_temporal->type() == "Constant");
}


/***********************************************************************//**
 * @brief Return spatial model component
 *
 * @return Pointer to spatial model component.
 *
 * Returns a pointer to the spatial model component of the model. The pointer
 * is of type GCTAModelSpatial. Note that a NULL pointer may be returned if
 * the model has no spatial model component.
 ***************************************************************************/
inline
GCTAModelSpatial* GCTAModelBackground::spatial(void) const
{
    return (m_spatial);
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
GModelSpectral* GCTAModelBackground::spectral(void) const
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
GModelTemporal* GCTAModelBackground::temporal(void) const
{
    return (m_temporal);
}

#endif /* GCTAMODELBACKGROUND_HPP */
