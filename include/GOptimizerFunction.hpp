/***************************************************************************
 *  GOptimizerFunction.hpp  -  Abstract base class for optimizer function  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GOptimizerFunction.hpp
 * @brief GOptimizerFunction abstract base class interface definition.
 * @author J. Knodlseder
 */

#ifndef GOPTIMIZERFUNCTION_HPP
#define GOPTIMIZERFUNCTION_HPP

/* __ Includes ___________________________________________________________ */
#include "GVector.hpp"


/***********************************************************************//**
 * @class GOptimizerFunction
 *
 * @brief GOptimizerFunction abstract base class interface defintion.
 ***************************************************************************/
class GOptimizerFunction {

public:
    // Iterator
    class iterator {
    friend class GOptimizerFunction;
    public:
        iterator();
        iterator(GOptimizerFunction *fct);
        ~iterator() { return; }
        iterator& operator++(void);          //!< Prefix
        iterator  operator++(int junk);      //!< Postfix
        bool      operator==(const iterator& it) const { return (m_index == it.m_index); }
        bool      operator!=(const iterator& it) const { return (m_index != it.m_index); }
        double    data(void) const { return m_base->get_data(); }
        double    model(void) const { return m_base->get_model(); }
        GVector   gradient(void) const { return m_base->get_gradient(); }
    protected:
        GOptimizerFunction* m_base;          //!< Pointer to base class
        int                 m_index;         //!< Iteration index
    };

    // Constructors and destructors
    GOptimizerFunction();
    GOptimizerFunction(const GOptimizerFunction& fct);
    virtual ~GOptimizerFunction();

    // Operators
    virtual GOptimizerFunction& operator= (const GOptimizerFunction& fct);

    // Virtual methods
    virtual void    first_item(void)   = 0;  //!< Move to first item
    virtual bool    next_item(void)    = 0;  //!< Move to next element (true if end)
    virtual double  get_data(void)     = 0;  //!< Get data value
    virtual double  get_model(void)    = 0;  //!< Get model value
    virtual GVector get_gradient(void) = 0;  //!< Get gradient
 
    // Implement methods
    iterator begin(void);
    iterator end(void);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GOptimizerFunction& fct);
    void free_members(void);

};

#endif /* GOPTIMIZERFUNCTION_HPP */
