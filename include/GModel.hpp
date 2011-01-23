/***************************************************************************
 *              GModel.hpp - Abstract virtual model base class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GModel.hpp
 * @brief GModel virtual base class interface definition.
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

/* __ Forward declarations _______________________________________________ */
class GEvent;
class GObservation;


/***********************************************************************//**
 * @class GModel
 *
 * @brief Model class interface defintion
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
    virtual GModelPar&       operator() (int index) = 0;
    virtual const GModelPar& operator() (int index) const = 0;
    GModel&                  operator= (const GModel& model);

    // Pure virtual methods
    virtual void        clear(void) = 0;
    virtual GModel*     clone(void) const = 0;
    virtual int         size(void) const = 0;
    virtual std::string type(void) const = 0;
    virtual double      eval(const GEvent& event, const GObservation& obs) = 0;
    virtual double      eval_gradients(const GEvent& event, const GObservation& obs) = 0;
    virtual void        read(const GXmlElement& xml) = 0;
    virtual void        write(GXmlElement& xml) const = 0;
    virtual std::string print(void) const = 0;

    // Implemented methods
    std::string name(void) const { return m_name; }
    void        name(const std::string& name) { m_name=name; }
    void        instruments(const std::string& instruments);
    std::string instruments(void) const;
    bool        isvalid(const std::string& name) const;

protected:
    // Protected methods
    void         init_members(void);
    void         copy_members(const GModel& model);
    void         free_members(void);
    std::string  print_name_instrument(void) const;
    virtual void set_pointers(void) = 0;

    // Proteced data members
    std::string              m_name;         //!< Model name
    std::vector<std::string> m_instruments;  //!< Instruments to which model applies
};

#endif /* GMODEL_HPP */
