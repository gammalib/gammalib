/***************************************************************************
 *      GCTAInstBackground.hpp - CTA instrument background model class     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014 by Juergen Knoedlseder                              *
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
 * @file GCTAInstBackground.hpp
 * @brief CTA instrument background model class definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTAINSTBACKGROUND_HPP
#define GCTAINSTBACKGROUND_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelData.hpp"
#include "GModelSpectral.hpp"
#include "GModelTemporal.hpp"
#include "GCTAEventList.hpp"

/* __ Forward declarations _______________________________________________ */


/***********************************************************************//**
 * @class GCTAInstBackground
 *
 * @brief CTA instrument background model class
 ***************************************************************************/
class GCTAInstBackground : public GModelData {

public:
    // Constructors and destructors
    GCTAInstBackground(void);
    GCTAInstBackground(const GCTAInstBackground& bgd);
    virtual ~GCTAInstBackground(void);

    // Operators
    GCTAInstBackground& operator=(const GCTAInstBackground& bgd);

    // Implemented pure virtual methods
    virtual void                clear(void);
    virtual GCTAInstBackground* clone(void) const;
    virtual std::string         type(void) const;
    virtual bool                is_constant(void) const;
    virtual double              eval(const GEvent& event,
                                     const GObservation& obs) const;
    virtual double              eval_gradients(const GEvent& event,
                                               const GObservation& obs) const;
    virtual double              npred(const GEnergy& obsEng,
                                      const GTime& obsTime,
                                      const GObservation& obs) const;
    virtual GCTAEventList*      mc(const GObservation& obs, GRan& ran) const;
    virtual void                read(const GXmlElement& xml);
    virtual void                write(GXmlElement& xml) const;
    virtual std::string         print(const GChatter& chatter = NORMAL) const;

    // Other methods
    GModelSpectral* spectral(void) const;
    GModelTemporal* temporal(void) const;
    
private:
    // Methods
    void            init_members(void);
    void            copy_members(const GCTAInstBackground& bgd);
    void            free_members(void);
    void            set_pointers(void);
    bool            valid_model(void) const;
    GModelSpectral* xml_spectral(const GXmlElement& spectral) const;
    GModelTemporal* xml_temporal(const GXmlElement& temporal) const;

    // Members
    GModelSpectral*   m_spectral;    //!< Spectral model
    GModelTemporal*   m_temporal;    //!< Temporal model
};


/***********************************************************************//**
 * @brief Return data model type
 *
 * @return Data model type.
 *
 * Returns the type of the data model.
 ***************************************************************************/
inline
std::string GCTAInstBackground::type(void) const
{
    return ("GCTAInstBackground");
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
bool GCTAInstBackground::is_constant(void) const
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
GModelSpectral* GCTAInstBackground::spectral(void) const
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
GModelTemporal* GCTAInstBackground::temporal(void) const
{
    return (m_temporal);
}

#endif /* GCTAINSTBACKGROUND_HPP */
