/***************************************************************************
 *               GCTAModelSkyCube.hpp - CTA sky cube model class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2020-2021 by Juergen Knoedlseder                         *
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
 * @file GCTAModelSkyCube.hpp
 * @brief CTA sky cube model class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTAMODELSKYCUBE_HPP
#define GCTAMODELSKYCUBE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GModelData.hpp"
#include "GFilename.hpp"
#include "GSkyMap.hpp"
#include "GEbounds.hpp"
#include "GNodeArray.hpp"
#include "GModelPar.hpp"
#include "GModelSpectral.hpp"
#include "GModelTemporal.hpp"
#include "GCTAEventList.hpp"

/* __ Forward declarations _______________________________________________ */
class GEvent;
class GObservation;
class GXmlElement;


/***********************************************************************//**
 * @class GCTAModelSkyCube
 *
 * @brief CTA sky cube model class
 *
 * This class implements a sky cube model for CTA.
 ***************************************************************************/
class GCTAModelSkyCube : public GModelData {

public:
    // Constructors and destructors
    GCTAModelSkyCube(void);
    explicit GCTAModelSkyCube(const GXmlElement& xml);
    GCTAModelSkyCube(const GFilename&      filename,
                     const GModelSpectral& spectral);
    GCTAModelSkyCube(const GFilename&      filename,
                     const GModelSpectral& spectral,
                     const GModelTemporal& temporal);
    GCTAModelSkyCube(const GCTAModelSkyCube& bgd);
    virtual ~GCTAModelSkyCube(void);

    // Operators
    virtual GCTAModelSkyCube& operator=(const GCTAModelSkyCube& model);

    // Implemented pure virtual methods
    virtual void              clear(void);
    virtual GCTAModelSkyCube* clone(void) const;
    virtual std::string       classname(void) const;
    virtual std::string       type(void) const;
    virtual bool              is_constant(void) const;
    virtual double            eval(const GEvent&       event,
                                   const GObservation& obs,
                                   const bool&         gradients = false) const;
    virtual double            npred(const GEnergy&       obsEng,
                                    const GTime&         obsTime,
                                    const GPolarization& obsPol,
                                    const GObservation&  obs) const;
    virtual GCTAEventList*    mc(const GObservation& obs, GRan& ran) const;
    virtual void              read(const GXmlElement& xml);
    virtual void              write(GXmlElement& xml) const;
    virtual std::string       print(const GChatter& chatter = NORMAL) const;
 
    // Other methods
    void             load(const GFilename& filename);
    const GFilename& filename(void) const;
    GModelSpectral*  spectral(void) const;
    GModelTemporal*  temporal(void) const;
    void             spectral(const GModelSpectral* spectral);
    void             temporal(const GModelTemporal* temporal);

protected:
    // Methods
    void            init_members(void);
    void            copy_members(const GCTAModelSkyCube& bgd);
    void            free_members(void);
    void            set_pointers(void);
    bool            valid_model(void) const;
    void            read_xml_spatial(const GXmlElement& xml);
    void            write_xml_spatial(GXmlElement& xml) const;
    GModelSpectral* xml_spectral(const GXmlElement& spectral) const;
    GModelTemporal* xml_temporal(const GXmlElement& temporal) const;

    // Members
    GFilename       m_filename;   //!< Filename
    GSkyMap         m_cube;       //!< Sky cube
    GEbounds        m_ebounds;    //!< Energy boundaries of sky cube
    GNodeArray      m_elogmeans;  //!< Mean log10(TeV) energies for the sky cube
    GModelPar       m_norm;       //!< Normalisation parameter
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
 * @return String containing the class name ("GCTAModelSkyCube").
 ***************************************************************************/
inline
std::string GCTAModelSkyCube::classname(void) const
{
    return ("GCTAModelSkyCube");
}


/***********************************************************************//**
 * @brief Return data model type
 *
 * @return Data model type.
 *
 * Returns the type of the data model.
 ***************************************************************************/
inline
std::string GCTAModelSkyCube::type(void) const
{
    return ("CTASkyCube");
}


/***********************************************************************//**
 * @brief Signals if sky model is temporally constant
 *
 * @return True if sky model is temporally constant, false otherwise.
 *
 * Signals if the sky model is temporally constant. A temporally constant
 * model is a model that has a temporal component of type "Constant". If
 * the model has no temporal component it is also consider being constant.
 ***************************************************************************/
inline
bool GCTAModelSkyCube::is_constant(void) const
{
    return ((m_temporal != NULL && m_temporal->type() == "Constant") ||
            (m_temporal == NULL));
}


/***********************************************************************//**
 * @brief Return filename
 *
 * @return Filename.
 ***************************************************************************/
inline
const GFilename& GCTAModelSkyCube::filename(void) const
{
    return (m_filename);
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
GModelSpectral* GCTAModelSkyCube::spectral(void) const
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
GModelTemporal* GCTAModelSkyCube::temporal(void) const
{
    return (m_temporal);
}

#endif /* GCTAMODELSKYCUBE_HPP */
