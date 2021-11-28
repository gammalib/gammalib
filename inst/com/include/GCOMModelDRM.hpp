/***************************************************************************
 *                GCOMModelDRM.hpp - COMPTEL DRM model class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2021 by Juergen Knoedlseder                              *
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
 * @file GCOMModelDRM.hpp
 * @brief COMPTEL DRM model class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GCOMMODELDRM_HPP
#define GCOMMODELDRM_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelData.hpp"
#include "GModelPar.hpp"
#include "GEvent.hpp"
#include "GObservation.hpp"
#include "GNodeArray.hpp"
#include "GXmlElement.hpp"
#include "GCOMEventCube.hpp"


/***********************************************************************//**
 * @class GCOMModelDRM
 *
 * @brief COMPTEL DRM model class
 *
 * This class implements a COMPTEL data space model. The model has a
 * normalization parameters.
 ***************************************************************************/
class GCOMModelDRM : public GModelData {

public:
    // Constructors and destructors
    GCOMModelDRM(void);
    explicit GCOMModelDRM(const GXmlElement& xml);
    GCOMModelDRM(const GCOMModelDRM& model);
    virtual ~GCOMModelDRM(void);

    // Operators
    virtual GCOMModelDRM& operator=(const GCOMModelDRM& model);

    // Implemented pure virtual methods
    virtual void           clear(void);
    virtual GCOMModelDRM*  clone(void) const;
    virtual std::string    classname(void) const;
    virtual std::string    type(void) const;
    virtual bool           is_constant(void) const;
    virtual double         eval(const GEvent& event,
                                const GObservation& obs,
                                const bool& gradients = false) const;
    virtual double         npred(const GEnergy& obsEng, const GTime& obsTime,
                                 const GObservation& obs) const;
    virtual GCOMEventCube* mc(const GObservation& obs, GRan& ran) const;
    virtual void           read(const GXmlElement& xml);
    virtual void           write(GXmlElement& xml) const;
    virtual std::string    print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCOMModelDRM& model);
    void free_members(void);

    // Proteced data members
    GModelPar m_norm;     //!< Normalisation
    GCOMDri   m_drm;      //!< Model
    GFilename m_filename; //!< Name of DRM FITS file
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCOMModelDRM").
 ***************************************************************************/
inline
std::string GCOMModelDRM::classname(void) const
{
    return ("GCOMModelDRM");
}


/***********************************************************************//**
 * @brief Return model type
 *
 * @return Model type.
 *
 * Returns the type of the model. The type for a DRM model is "DRM".
 ***************************************************************************/
inline
std::string GCOMModelDRM::type(void) const
{
    return ("DRM");
}


/***********************************************************************//**
 * @brief Signals if model is temporally constant
 *
 * @return True if model is temporally constant, false otherwise.
 *
 * Signals if the model is temporally constant. By definition, a DRM model
 * is always temporally constant.
 ***************************************************************************/
inline
bool GCOMModelDRM::is_constant(void) const
{
    return (true);
}

#endif /* GCOMMODELDRM_HPP */
