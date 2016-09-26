/***************************************************************************
 *      GCTAModelCubeBackground.hpp - CTA cube background model class      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013-2016 by Michael Mayer                               *
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
 * @file GCTAModelCubeBackground.hpp
 * @brief CTA cube background model class interface definition
 * @author Michael Mayer
 */

#ifndef GCTAMODELCUBEBACKGROUND_HPP
#define GCTAMODELCUBEBACKGROUND_HPP

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
#include "GCTAObservation.hpp"
#include "GCTAEventList.hpp"
#include "GModelSpatial.hpp"


/***********************************************************************//**
 * @class GCTAModelCubeBackground
 *
 * @brief CTA cube background model class
 *
 * This class implements a cube background model for CTA.
 ***************************************************************************/
class GCTAModelCubeBackground : public GModelData {

public:
    // Constructors and destructors
    GCTAModelCubeBackground(void);
    explicit GCTAModelCubeBackground(const GXmlElement& xml);
    explicit GCTAModelCubeBackground(const GModelSpectral& spectral);
    GCTAModelCubeBackground(const GCTAModelCubeBackground& bgd);
    virtual ~GCTAModelCubeBackground(void);

    // Operators
    virtual GCTAModelCubeBackground& operator=(const GCTAModelCubeBackground& model);

    // Implemented pure virtual methods
    virtual void                     clear(void);
    virtual GCTAModelCubeBackground* clone(void) const;
    virtual std::string              classname(void) const;
    virtual std::string              type(void) const;
    virtual bool                     is_constant(void) const;
    virtual double                   eval(const GEvent& event,
                                          const GObservation& obs,
                                          const bool& gradients = false) const;
    virtual double                   npred(const GEnergy& obsEng,
                                           const GTime& obsTime,
                                           const GObservation& obs) const;
    virtual GCTAEventList*           mc(const GObservation& obs, GRan& ran) const;
    virtual void                     read(const GXmlElement& xml);
    virtual void                     write(GXmlElement& xml) const;
    virtual std::string              print(const GChatter& chatter = NORMAL) const; 
 
    // Other methods
    GModelSpectral* spectral(void) const;
    GModelTemporal* temporal(void) const;

protected:
    // Methods
    void            init_members(void);
    void            copy_members(const GCTAModelCubeBackground& bgd);
    void            free_members(void);
    void            set_pointers(void);
    bool            valid_model(void) const;
    GModelSpectral* xml_spectral(const GXmlElement& spectral) const;
    GModelTemporal* xml_temporal(const GXmlElement& temporal) const;

    // Members
    GModelSpectral* m_spectral;   //!< Spectral model
    GModelTemporal* m_temporal;   //!< Temporal model

    // Npred cache
    mutable std::vector<std::string> m_npred_names;    //!< Model names
    mutable std::vector<GEnergy>     m_npred_energies; //!< Model energy
    mutable std::vector<GTime>       m_npred_times;    //!< Model time
    mutable std::vector<double>      m_npred_values;   //!< Model values
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCTAModelCubeBackground").
 ***************************************************************************/
inline
std::string GCTAModelCubeBackground::classname(void) const
{
    return ("GCTAModelCubeBackground");
}


/***********************************************************************//**
 * @brief Return data model type
 *
 * @return Data model type.
 *
 * Returns the type of the data model.
 ***************************************************************************/
inline
std::string GCTAModelCubeBackground::type(void) const
{
    return ("CTACubeBackground");
}


/***********************************************************************//**
 * @brief Signals if sky model is temporally constant
 *
 * @return True if sky model is temporally constant, false otherwise.
 *
 * Signals if the sky model is temporally constant. A temporally constant
 * model is a model that has a temporal component of type "Constant".
 ***************************************************************************/
inline
bool GCTAModelCubeBackground::is_constant(void) const
{
    return (m_temporal != NULL && m_temporal->type() == "Constant");
}


/***********************************************************************//**
 * @brief Return spectral model component
 *
 * @return Pointer to spectral model component.
 *
 * Returns a pointer to the spectral model component of the model. The
 * pointer is of type GModelSpectral. Note that a NULL pointer may be
 * returned if the sky model has no spectral model component.
 ***************************************************************************/
inline
GModelSpectral* GCTAModelCubeBackground::spectral(void) const
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
 * returned if the sky model has no temporal model component.
 ***************************************************************************/
inline
GModelTemporal* GCTAModelCubeBackground::temporal(void) const
{
    return (m_temporal);
}

#endif /* GCTAMODELCUBEBACKGROUND_HPP */
