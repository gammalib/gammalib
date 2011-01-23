/***************************************************************************
 *         GModel.i - Abstract virtual model base class python I/F         *
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
 * @file GModel.i
 * @brief GModel class python interface
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GModel.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GModel
 *
 * @brief GModel class python interface defintion
 ***************************************************************************/
class GModel {
public:
    // Constructors and destructors
    GModel(void);
    explicit GModel(const GXmlElement& xml);
    GModel(const GModel& model);
    virtual ~GModel(void);

    // Pure virtual methods
    virtual void        clear(void) = 0;
    virtual GModel*     clone(void) const = 0;
    virtual int         size(void) const = 0;
    virtual std::string type(void) const = 0;
    virtual double      eval(const GEvent& event, const GObservation& obs) = 0;
    virtual double      eval_gradients(const GEvent& event, const GObservation& obs) = 0;
    virtual void        read(const GXmlElement& xml) = 0;
    virtual void        write(GXmlElement& xml) const = 0;

    // Implemented methods
    std::string name(void) const;
    void        name(const std::string& name);
    void        instruments(const std::string& instruments);
    std::string instruments(void) const;
    bool        isvalid(const std::string& name) const;

protected:
    virtual void set_pointers(void) = 0;
};


/***********************************************************************//**
 * @brief GModel class extension
 ***************************************************************************/
%extend GModel {
    char *__str__() {
        return tochar(self->print());
    }
    GModelPar __getitem__(int index) {
    if (index >= 0 && index < self->size())
        return (*self)(index);
    else
        throw GException::out_of_range("__getitem__(int)", index, self->size());
    }
    void __setitem__(int index, const GModelPar& val) {
        if (index>=0 && index < self->size())
            (*self)(index) = val;
        else
            throw GException::out_of_range("__setitem__(int)", index, self->size());
    }
};
