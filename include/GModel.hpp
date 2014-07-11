/***************************************************************************
 *              GModel.hpp - Abstract virtual model base class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2014 by Juergen Knoedlseder                         *
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
#include "GXmlElement.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"

/* __ Forward declarations _______________________________________________ */
class GEvent;
class GObservation;


/***********************************************************************//**
 * @class GModel
 *
 * @brief Abstract model class
 *
 * This class implements a parametric model. The class has two methods for
 * model evaluation: eval() and eval_gradients().
 * The eval() method evaluates the model for a given event and observation.
 * In addition, eval_gradients() also sets the parameter gradients of the
 * model.
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
    virtual std::string type(void) const = 0;
    virtual bool        is_constant(void) const = 0;
    virtual double      eval(const GEvent& event,
                             const GObservation& obs) const = 0;
    virtual double      eval_gradients(const GEvent& event,
                                       const GObservation& obs) const = 0;
    virtual double      npred(const GEnergy& obsEng, const GTime& obsTime,
                              const GObservation& obs) const = 0;
    virtual void        read(const GXmlElement& xml) = 0;
    virtual void        write(GXmlElement& xml) const = 0;
    virtual std::string print(const GChatter& chatter = NORMAL) const = 0;

    // Implemented methods
    int                 size(void) const;
    GModelPar&          at(const int& index);
    const GModelPar&    at(const int& index) const;
    bool                has_par(const std::string& name) const;
    const std::string&  name(void) const;
    void                name(const std::string& name);
    const double&       ts(void) const;
    void                ts(const double& ts);
    std::string         instruments(void) const;
    void                instruments(const std::string& instruments);
    GModelPar           scale(const std::string& instrument) const;
    void                scale(const GModelPar& par);
    std::string         ids(void) const;
    void                ids(const std::string& ids);
    bool                is_valid(const std::string& instruments,
                                const std::string& ids) const;

protected:
    // Protected methods
    void         init_members(void);
    void         copy_members(const GModel& model);
    void         free_members(void);
    void         read_scales(const GXmlElement& xml);
    void         write_scales(GXmlElement& xml) const;
    std::string  print_attributes(void) const;

    // Protected members
    std::string              m_name;         //!< Model name
    std::vector<std::string> m_instruments;  //!< Instruments to which model applies
    std::vector<GModelPar>   m_scales;       //!< Model instrument scale factors
    std::vector<std::string> m_ids;          //!< Identifiers to which model applies
    std::vector<GModelPar*>  m_pars;         //!< Pointers to all model parameters
    double                   m_ts;           //!< Test Statistic of the Model
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
    return (m_pars.size());
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
 * @brief Return TS value
 *
 * @return TS value.
 *
 * Returns the TS value.
 ***************************************************************************/
inline
const double& GModel::ts(void) const
{
    return m_ts;
}


/***********************************************************************//**
 * @brief Set TS value
 *
 * @param[in] TS value.
 *
 * Set the TS value.
 ***************************************************************************/
inline
void GModel::ts(const double& ts)
{
    m_ts = ts;
    return;
}


#endif /* GMODEL_HPP */
