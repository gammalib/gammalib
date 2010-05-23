/***************************************************************************
 *             GObservations.hpp  -  Observation container class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2010 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GObservations.hpp
 * @brief GObservations container class interface definition.
 * @author J. Knodlseder
 */

#ifndef GOBSERVATIONS_HPP
#define GOBSERVATIONS_HPP

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include "GObservation.hpp"
#include "GEvent.hpp"
#include "GOptimizer.hpp"
#include "GOptimizerFunction.hpp"
#include "GModels.hpp"


/***********************************************************************//**
 * @class GObservations
 *
 * @brief GObservations container class interface defintion.
 ***************************************************************************/
class GObservations {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GObservations& obs);

public:
    // Constructors and destructors
    GObservations(void);
    GObservations(const GObservations& obs);
    virtual ~GObservations(void);

    // Operators
    GObservations& operator= (const GObservations& obs);

    // Methods
	void          append(GObservation &obs);
    void          models(const GModels& models) { m_models=models; return; }
    int           elements(void) const { return m_num; }
    GObservation* observation(int index) const;
    GModels*      models(void) { return &m_models; }
    void          optimize(GOptimizer& opt);

    // Event iterator
    class iterator {
    friend class GObservations;
    public:
        iterator();
        iterator(GObservations *obs);
        ~iterator() { return; }
        iterator& operator++(void);                // Prefix
        iterator  operator++(int junk);            // Postfix
        bool      operator==(const iterator& it) const 
                  { return ((m_index == it.m_index) && (m_event == it.m_event)); }
        bool      operator!=(const iterator& it) const
                  { return ((m_index != it.m_index) || (m_event != it.m_event)); }
        GEvent&   operator*(void) { return *m_event; }
        GEvent*   operator->(void) { return &(*m_event); }
    protected:
        int               m_index;    //!< Actual observation index [0,m_num-1]
        GEvents::iterator m_event;    //!< Iterator on actual event
        GEvents::iterator m_end;      //!< Iterator on observation end
        GObservation*     m_obs;      //!< Pointer to actual observation
        GObservations*    m_data;     //!< Pointer to GObservations object
    };
    iterator begin(void);
    iterator end(void);
    
    // Optimizer
    class optimizer : public GOptimizerFunction {
    public:
        optimizer();
        optimizer(GObservations *obs);
        optimizer(const optimizer& fct);
        ~optimizer();
        optimizer& operator= (const optimizer& fct);
        void           eval(const GOptimizerPars& pars);
        double*        value(void) { return &m_value; }
        GVector*       gradient(void) { return m_gradient; }
        GSparseMatrix* covar(void) { return m_covar; }
    protected:
        void           init_members(void);
        void           copy_members(const optimizer& fct);
        void           free_members(void);
        double         m_value;       //!< Function value
        GVector*       m_gradient;    //!< Pointer to gradient vector
        GSparseMatrix* m_covar;       //!< Pointer to covariance matrix
        GObservations* m_data;        //!< Pointer to GObservations object
    };

protected:
    // Protected methods
    void           init_members(void);
    void           copy_members(const GObservations& obs);
    void           free_members(void);

    // Protected data area
	int            m_num;            //!< Number of observations
	GObservation** m_obs;            //!< Pointers to observations
    GModels        m_models;         //!< Models
	
};

#endif /* GOBSERVATIONS_HPP */
