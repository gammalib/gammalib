/***************************************************************************
 *              GModel.hpp - Abstract virtual model base class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2020 by Juergen Knoedlseder                         *
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
 * @file GModel.hpp
 * @brief Abstract model base class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GMODEL_HPP
#define GMODEL_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include "GBase.hpp"
#include "GModelPar.hpp"
#include "GModelAssociations.hpp"
#include "GXmlElement.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"

/* __ Forward declarations _______________________________________________ */
class GVector;
class GMatrixSparse;
class GEvent;
class GObservation;


/***********************************************************************//**
 * @class GModel
 *
 * @brief Abstract model class
 *
 * This class implements a parametric model. The eval() method evaluates the
 * model for a given event and observation. If the gradients parameter is set
 * to true, the eval() method also computes analytical parameter gradients if
 * they exist.
 *
 * A model has the following attributes:
 * - @p name
 * - @p type
 * - @p instruments
 * - @p ids
 * - @p scales
 *
 * The model @p name is a text string that names the model. The name()
 * methods allow setting and retrieving the model name. The model handling
 * does not actual depend on the value of this text string.
 *
 * The model @p type is a text string that specifies the kind of the model.
 * Examples are @p PointSource, @p ExtendedSource or @p DiffuseSource for
 * sky models. The type() method allows retrieving the model type. The model
 * handling does not actual depend on the value of this text string.
 *
 * The model @p instruments is a list of text strings that specifies the
 * instruments to which the model applies. The is_valid() method will check
 * a given instrument name against that list to verify if the model applies
 * to an instrument. The instruments() method returns a comma separated
 * list of all instruments. Another instance of this method allows setting
 * the instrument list from a comma separated list. If the @p instruments
 * list is empty, the model applies to all instruments. 
 *
 * The model @p ids is a list of text strings that specifies the observation
 * identifiers to which the model applies. This allows specifying of models
 * for a specific list of observations. The is_valid() method will check
 * a given observation identifier against that list to verify if the model
 * applies. The ids() method returns a comma separated list of all
 * observation identifiers. Another instance of this method allows setting
 * the observation identifiers from a comma separated list. If the @p ids
 * list is empty, the model applies to all observation identifiers.
 *
 * The model @p scales is a list of model parameters that specify an
 * instrument dependent scaling factor. This allows to prescale a model for
 * a given instrument. The scale() methods allow to set and to retrieve the
 * scale factor for a specific instrument.
 *
 * As the class holds simply a collection of model parameters, it should
 * neither deal with allocation and deallocation, nor with cloning of
 * model parameters. This will be done by the classes that actually
 * implement the model parameters.
 ***************************************************************************/
class GModel : public GBase {

public:
    // Constructors and destructors
    GModel(void);
    explicit GModel(const GXmlElement& xml);
    GModel(const GModel& model);
    virtual ~GModel(void);

    // Operators
    virtual GModel&          operator=(const GModel& model);
    virtual GModelPar&       operator[](const int& index);
    virtual const GModelPar& operator[](const int& index) const;
    virtual GModelPar&       operator[](const std::string& name);
    virtual const GModelPar& operator[](const std::string& name) const;

    // Pure virtual methods
    virtual void        clear(void) = 0;
    virtual GModel*     clone(void) const = 0;
    virtual std::string classname(void) const = 0;
    virtual std::string type(void) const = 0;
    virtual bool        is_constant(void) const = 0;
    virtual double      eval(const GEvent& event,
                             const GObservation& obs,
                             const bool& gradients = false) const = 0;
    virtual GVector     eval(const GObservation& obs,
                             GMatrixSparse* gradients = NULL) const = 0;
    virtual double      npred(const GEnergy& obsEng, const GTime& obsTime,
                              const GObservation& obs) const = 0;
    virtual void        read(const GXmlElement& xml) = 0;
    virtual void        write(GXmlElement& xml) const = 0;
    virtual std::string print(const GChatter& chatter = NORMAL) const = 0;

    // Implemented methods
    int                       size(void) const;
    int                       scales(void) const;
    GModelPar&                at(const int& index);
    const GModelPar&          at(const int& index) const;
    bool                      has_par(const std::string& name) const;
    bool                      has_scales(void) const;
    const std::string&        name(void) const;
    void                      name(const std::string& name);
    const double&             ts(void) const;
    void                      ts(const double& ts);
    const bool&               tscalc(void) const;
    void                      tscalc(const bool& tscalc);
    const bool&               has_ts(void) const;
    std::string               instruments(void) const;
    void                      instruments(const std::string& instruments);
    GModelPar&                scale(const int& index);
    const GModelPar&          scale(const int& index) const;
    GModelPar                 scale(const std::string& instrument) const;
    void                      scale(const GModelPar& par);
    std::string               ids(void) const;
    void                      ids(const std::string& ids);
    bool                      is_valid(const std::string& instruments,
                                       const std::string& ids) const;
    const GModelAssociations& associations(void) const;
    void                      associations(const GModelAssociations& associations);
    const bool&               has_eval_indices(void) const;
    const std::vector<int>&   eval_indices(void) const;

protected:
    // Protected methods
    void         init_members(void);
    void         copy_members(const GModel& model);
    void         free_members(void);
    void         read_attributes(const GXmlElement& xml);
    void         write_attributes(GXmlElement& xml) const;
    std::string  print_attributes(void) const;
    void         read_scales(const GXmlElement& xml);
    void         write_scales(GXmlElement& xml) const;

    // Protected members
    std::string              m_name;         //!< Model name
    std::vector<std::string> m_instruments;  //!< Instruments to which model applies
    std::vector<GModelPar>   m_scales;       //!< Model instrument scale factors
    std::vector<std::string> m_ids;          //!< Identifiers to which model applies
    std::vector<GModelPar*>  m_pars;         //!< Pointers to all model parameters
    GModelAssociations       m_associations; //!< Model associations
    bool                     m_has_ts;       //!< Signals if TS is available
    bool                     m_has_tscalc;   //!< Signals if tscalc attribute is available
    bool                     m_tscalc;       //!< Signals if TS should be computed
    double                   m_ts;           //!< Test Statistic of the model

    // Protected members for efficient gradient access
    mutable bool             m_has_eval_inx; //!< Signal that parameter indices
                                             //!< updated by eval() method exist
    mutable std::vector<int> m_eval_inx;     //!< List of parameter indices updated
                                             //!< by eval() method
};


/***********************************************************************//**
 * @brief Returns reference to model parameter by index
 *
 * @param[in] index Parameter index [0,...,size()-1].
 * @return Reference to model parameter.
 *
 * Returns a reference to the model parameter of the specified @p index.
 ***************************************************************************/
inline
GModelPar& GModel::operator[](const int& index)
{
    // Return reference
    return *(m_pars[index]);
}


/***********************************************************************//**
 * @brief Returns reference to model parameter by index (const version)
 *
 * @param[in] index Parameter index [0,...,size()-1].
 * @return Const reference to model parameter.
 *
 * Returns a const reference to the model parameter of the specified
 ***************************************************************************/
inline
const GModelPar& GModel::operator[](const int& index) const
{
    // Return reference
    return *(m_pars[index]);
}


/***********************************************************************//**
 * @brief Return number of parameters in model
 *
 * @return Number of parameters in model.
 *
 * Returns the number of parameters in the model.
 ***************************************************************************/
inline
int GModel::size(void) const
{
    return (int)m_pars.size();
}


/***********************************************************************//**
 * @brief Return number of scale parameters in model
 *
 * @return Number of scale parameters in model.
 *
 * Returns the number of scale parameters in the model.
 ***************************************************************************/
inline
int GModel::scales(void) const
{
    return (int)m_scales.size();
}


/***********************************************************************//**
 * @brief Return parameter name
 *
 * @return Parameter name.
 *
 * Returns the parameter name.
 ***************************************************************************/
inline
const std::string& GModel::name(void) const
{
    return m_name;
}


/***********************************************************************//**
 * @brief Set parameter name
 *
 * @param[in] name Parameter name.
 *
 * Set the parameter name.
 ***************************************************************************/
inline
void GModel::name(const std::string& name)
{
    m_name = name;
    return;
}


/***********************************************************************//**
 * @brief Return Test Statistic value
 *
 * @return Test Statistic value.
 *
 * Returns the Test Statistic value. The Test Statistic value is twice the
 * difference between the log-likelihood of the null hypothesis and the
 * alternative hypothesis.
 ***************************************************************************/
inline
const double& GModel::ts(void) const
{
    return m_ts;
}


/***********************************************************************//**
 * @brief Set Test Statistic value
 *
 * @param[in] ts Test Statistic value.
 *
 * Set the Test Statistic value. The Test Statistic value is twice the
 * difference between the log-likelihood of the null hypothesis and the
 * alternative hypothesis.
 ***************************************************************************/
inline
void GModel::ts(const double& ts)
{
    m_ts     = ts;
    m_has_ts = true; //!< Signals that TS is now available
    return;
}


/***********************************************************************//**
 * @brief Return Test Statistic computation flag
 *
 * @return Test Statistic computation flag; true is Test Statistic should be
 *         computed.
 *
 * Returns the flag that signals whether Test Statistic values should be
 * computed.
 ***************************************************************************/
inline
const bool& GModel::tscalc(void) const
{
    return m_tscalc;
}


/***********************************************************************//**
 * @brief Set Test Statistic computation flag
 *
 * @param[in] tscalc Test Statistic computation flag.
 *
 * Set the flag that signals whether Test Statistic values should be
 * computed.
 ***************************************************************************/
inline
void GModel::tscalc(const bool& tscalc)
{
    m_tscalc     = tscalc;
    m_has_tscalc = true; //!< Signals that tscalc is now available
    return;
}


/***********************************************************************//**
 * @brief Signals that model has scales
 *
 * @return True if model has scale factors.
 ***************************************************************************/
inline
bool GModel::has_scales(void) const
{
    return (!m_scales.empty());
}


/***********************************************************************//**
 * @brief Signals that model has Test Statistics value
 *
 * @return True if model has TS value.
 ***************************************************************************/
inline
const bool& GModel::has_ts(void) const
{
    return m_has_ts;
}


/***********************************************************************//**
 * @brief Return model associations
 *
 * @return Model associations.
 *
 * Returns model associations.
 ***************************************************************************/
inline
const GModelAssociations& GModel::associations(void) const
{
    return m_associations;
}


/***********************************************************************//**
 * @brief Set model associations
 *
 * @param[in] associations Model associations.
 *
 * Set the model associations.
 ***************************************************************************/
inline
void GModel::associations(const GModelAssociations& associations)
{
    m_associations = associations;
    return;
}


/***********************************************************************//**
 * @brief Signals that parameter indices updated by eval() method were set
 *
 * @return True if parameter indices updated by eval() method were set.
 ***************************************************************************/
inline
const bool& GModel::has_eval_indices(void) const
{
    return m_has_eval_inx;
}


/***********************************************************************//**
 * @brief Return vector of parameter indices updated by eval() method
 *
 * @return Vector of parameter indices updated by eval() method.
 ***************************************************************************/
inline
const std::vector<int>& GModel::eval_indices(void) const
{
    return m_eval_inx;
}

#endif /* GMODEL_HPP */
