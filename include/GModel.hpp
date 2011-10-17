/***************************************************************************
 *              GModel.hpp - Abstract virtual model base class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2011 by Jurgen Knodlseder                           *
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
 * @author J. Knodlseder
 */

#ifndef GMODEL_HPP
#define GMODEL_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include <iostream>
#include "GLog.hpp"
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
 * This abstract virtual base class implements the methods that handle the
 * model name and the applicable instruments. It also implements friend
 * operators for logging.
 *
 * As the class holds simpliy a collection of model parameters, it should
 * neither deal with allocation and deallocation, nor with cloning of
 * model parameters. This will be done by the classes that actually
 * implement the model parameters.
 ***************************************************************************/
class GModel {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GModel& model);
    friend GLog&         operator<< (GLog& log, const GModel& model);

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
    virtual double      eval(const GEvent& event,
                             const GObservation& obs) const = 0;
    virtual double      eval_gradients(const GEvent& event,
                                       const GObservation& obs) const = 0;
    virtual double      npred(const GEnergy& obsEng, const GTime& obsTime,
                              const GObservation& obs) const = 0;
    virtual void        read(const GXmlElement& xml) = 0;
    virtual void        write(GXmlElement& xml) const = 0;
    virtual std::string print(void) const = 0;

    // Implemented methods
    int                 size(void) const { return m_pars.size(); }
    std::string         name(void) const { return m_name; }
    void                name(const std::string& name) { m_name=name; }
    void                instruments(const std::string& instruments);
    std::string         instruments(void) const;
    bool                isvalid(const std::string& name) const;

protected:
    // Protected methods
    void         init_members(void);
    void         copy_members(const GModel& model);
    void         free_members(void);
    std::string  print_name_instrument(void) const;

    // Proteced members
    std::string              m_name;         //!< Model name
    std::vector<std::string> m_instruments;  //!< Instruments to which model applies
    std::vector<GModelPar*>  m_pars;         //!< Pointers to all model parameters
};

#endif /* GMODEL_HPP */
