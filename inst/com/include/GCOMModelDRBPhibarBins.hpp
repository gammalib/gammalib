/***************************************************************************
 * GCOMModelDRBPhibarBins.hpp - COMPTEL DRB model Phibar bin fitting class *
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
 * @file GCOMModelDRBPhibarBins.hpp
 * @brief COMPTEL DRB model Phibar bin fitting class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GCOMMMODELDRBPHIBARBINS_HPP
#define GCOMMMODELDRBPHIBARBINS_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GModelData.hpp"
#include "GModelPar.hpp"
#include "GCOMEventCube.hpp"

/*__ Forward declarations _________________________________________________ */
class GEvent;
class GObservation;
class GXmlElement;


/***********************************************************************//**
 * @class GCOMModelDRBPhibarBins
 *
 * @brief COMPTEL DRB model Phibar bin fitting class
 *
 * This class implements a COMPTEL background model that is based on fitting
 * of a DRB model. The model parameters are a set of normalization parameters
 * defined for a given Phibar value.
 ***************************************************************************/
class GCOMModelDRBPhibarBins : public GModelData {

public:
    // Constructors and destructors
    GCOMModelDRBPhibarBins(void);
    explicit GCOMModelDRBPhibarBins(const GXmlElement& xml);
    GCOMModelDRBPhibarBins(const GCOMModelDRBPhibarBins& model);
    virtual ~GCOMModelDRBPhibarBins(void);

    // Operators
    virtual GCOMModelDRBPhibarBins& operator=(const GCOMModelDRBPhibarBins& model);

    // Implemented pure virtual methods
    virtual void                    clear(void);
    virtual GCOMModelDRBPhibarBins* clone(void) const;
    virtual std::string             classname(void) const;
    virtual std::string             type(void) const;
    virtual bool                    is_constant(void) const;
    virtual double                  eval(const GEvent& event,
                                         const GObservation& obs,
                                         const bool& gradients = false) const;
    virtual double                  npred(const GEnergy& obsEng, const GTime& obsTime,
                                          const GObservation& obs) const;
    virtual GCOMEventCube*          mc(const GObservation& obs, GRan& ran) const;
    virtual void                    read(const GXmlElement& xml);
    virtual void                    write(GXmlElement& xml) const;
    virtual std::string             print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCOMModelDRBPhibarBins& model);
    void free_members(void);
    void set_pointers(void);

    // Proteced data members
    std::vector<GModelPar> m_values; //!< Normalisation values
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCOMModelDRBPhibarBins").
 ***************************************************************************/
inline
std::string GCOMModelDRBPhibarBins::classname(void) const
{
    return ("GCOMModelDRBPhibarBins");
}


/***********************************************************************//**
 * @brief Return model type
 *
 * @return Model type.
 *
 * Returns the type of the model. The type for a DRB Phibar bin fitting model
 * is "DRBPhibarBins".
 ***************************************************************************/
inline
std::string GCOMModelDRBPhibarBins::type(void) const
{
    return ("DRBPhibarBins");
}


/***********************************************************************//**
 * @brief Signals if model is temporally constant
 *
 * @return True if model is temporally constant, false otherwise.
 *
 * Signals if the model is temporally constant. By definition, a DRB Phibar
 * bin fitting model is always temporally constant.
 ***************************************************************************/
inline
bool GCOMModelDRBPhibarBins::is_constant(void) const
{
    return (true);
}

#endif /* GCOMMMODELDRBPHIBARBINS_HPP */
