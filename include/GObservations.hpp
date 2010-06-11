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
 *
 * This is the main user interface class that provides high level access to
 * gamma-ray observations and their analysis. Observations are appended using
 * the append() method and can be access using the operator(). The class
 * implements an event iterator GObservations::iterator that allows iterating
 * over all events in all observations.
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
    GObservations&      operator= (const GObservations& obs);
    GObservation&       operator() (int index);
    const GObservation& operator() (int index) const;

    // Methods
    void     append(GObservation &obs);
    int      size(void) const { return m_num; }
    void     models(const GModels& models) { m_models=models; return; }
    GModels* models(void) { return &m_models; }
    void     optimize(GOptimizer& opt);

    // Event iterator
    class iterator {
    friend class GObservations;
    public:
        // Constructors and destructors
        iterator(void);
        iterator(GObservations *obs);
        ~iterator(void) { return; }

        // Operators
        iterator&     operator++(void);                // Prefix
        iterator      operator++(int junk);            // Postfix
        bool          operator==(const iterator& it) const 
                      { return ((m_index == it.m_index) && (m_event == it.m_event)); }
        bool          operator!=(const iterator& it) const
                      { return ((m_index != it.m_index) || (m_event != it.m_event)); }
        GEvent&       operator*(void) { return *m_event; }
        GEvent*       operator->(void) { return &(*m_event); }

        // Methods
        GObservation* obs(void) { return m_obs; }
    protected:
        int               m_index;    //!< Actual observation index [0,m_num-1]
        GEvents::iterator m_event;    //!< Iterator on actual event
        GEvents::iterator m_end;      //!< Iterator on observation end
        GObservation*     m_obs;      //!< Pointer to actual observation
        GObservations*    m_this;     //!< Pointer to GObservations object
    };
    iterator begin(void);
    iterator end(void);

    // Optimizer
    class optimizer : public GOptimizerFunction {
    public:
        // Constructors and destructors
        optimizer(void);
        optimizer(GObservations *obs);
        optimizer(const optimizer& fct);
        ~optimizer(void);

        // Operators
        optimizer& operator= (const optimizer& fct);

        // Methods
        void           eval(const GOptimizerPars& pars);
        double*        value(void) { return &m_value; }
        GVector*       gradient(void) { return m_gradient; }
        GSparseMatrix* covar(void) { return m_covar; }
    protected:
        void           init_members(void);
        void           copy_members(const optimizer& fct);
        void           free_members(void);
        double         m_value;       //!< Function value
        double         m_npred;       //!< Total number of predicted events
        double         m_minmod;      //!< Minimum model value
        GVector*       m_gradient;    //!< Pointer to gradient vector
        GSparseMatrix* m_covar;       //!< Pointer to covariance matrix
        GObservations* m_this;        //!< Pointer to GObservations object
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
