/***************************************************************************
 *       GCOMModelDRBPhibarNodes.hpp - COMPTEL DRB model fitting class     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2022 by Juergen Knoedlseder                         *
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
 * @file GCOMModelDRBPhibarNodes.hpp
 * @brief COMPTEL DRB Phibar nodes model fitting class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GCOMMODELDRBPHIBARNODES_HPP
#define GCOMMODELDRBPHIBARNODES_HPP

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
 * @class GCOMModelDRBPhibarNodes
 *
 * @brief COMPTEL DRB Phibar nodes model fitting class
 *
 * This class implements a COMPTEL background model that is based on fitting
 * of a DRB model. The model parameters are a set of normalization parameters
 * defined for a given Phibar value, and linear interpolations are performed
 * between these values.
 ***************************************************************************/
class GCOMModelDRBPhibarNodes : public GModelData {

public:
    // Constructors and destructors
    GCOMModelDRBPhibarNodes(void);
    explicit GCOMModelDRBPhibarNodes(const GXmlElement& xml);
    GCOMModelDRBPhibarNodes(const GCOMModelDRBPhibarNodes& model);
    virtual ~GCOMModelDRBPhibarNodes(void);

    // Operators
    virtual GCOMModelDRBPhibarNodes& operator=(const GCOMModelDRBPhibarNodes& model);

    // Implemented pure virtual methods
    virtual void                     clear(void);
    virtual GCOMModelDRBPhibarNodes* clone(void) const;
    virtual std::string              classname(void) const;
    virtual std::string              type(void) const;
    virtual bool                     is_constant(void) const;
    virtual double                   eval(const GEvent& event,
                                          const GObservation& obs,
                                          const bool& gradients = false) const;
    virtual double                   npred(const GEnergy& obsEng, const GTime& obsTime,
                                           const GObservation& obs) const;
    virtual GCOMEventCube*           mc(const GObservation& obs, GRan& ran) const;
    virtual void                     read(const GXmlElement& xml);
    virtual void                     write(GXmlElement& xml) const;
    virtual std::string              print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCOMModelDRBPhibarNodes& model);
    void free_members(void);
    void set_pointers(void);
    void set_cache(void) const;
    void update_cache(void) const;

    // Proteced data members
    std::vector<GModelPar>      m_phibars;      //!< Node Phibar values
    std::vector<GModelPar>      m_values;       //!< Node values

    // Evaluation cache
    mutable bool                m_scale;        //!< Model is a scale factor
    mutable bool                m_fixed;        //!< All Phibar values are fixed
    mutable std::vector<double> m_old_phibars;  //!< Old Phibar values
    mutable GNodeArray          m_nodes;        //!< Phibar node values
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCOMModelDRBPhibarNodes").
 ***************************************************************************/
inline
std::string GCOMModelDRBPhibarNodes::classname(void) const
{
    return ("GCOMModelDRBPhibarNodes");
}


/***********************************************************************//**
 * @brief Return model type
 *
 * @return Model type.
 *
 * Returns the type of the model. The type for a DRB fitting model is
 * "DRBPhibarNodes".
 ***************************************************************************/
inline
std::string GCOMModelDRBPhibarNodes::type(void) const
{
    return ("DRBPhibarNodes");
}


/***********************************************************************//**
 * @brief Signals if model is temporally constant
 *
 * @return True if model is temporally constant, false otherwise.
 *
 * Signals if the model is temporally constant. By definition, a DRB fitting
 * model is always temporally constant.
 ***************************************************************************/
inline
bool GCOMModelDRBPhibarNodes::is_constant(void) const
{
    return (true);
}

#endif /* GCOMMODELDRBPHIBARNODES_HPP */
