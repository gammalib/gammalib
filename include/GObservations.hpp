/***************************************************************************
 *             GObservations.hpp  -  Observation container class           *
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
 * @file GObservations.hpp
 * @brief Observations container class interface definition
 * @author J. Knodlseder
 */

#ifndef GOBSERVATIONS_HPP
#define GOBSERVATIONS_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include <iostream>
#include "GLog.hpp"
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
 *
 * @todo Implement models pointer instead of models elements since models
 * should be attached without loosing the external link to the model
 * parameters.
 *
 * @todo Implement append() method using a pointer instead of a reference
 * to make clear that the observation is linked from the outside (and thus
 * has to be preserved). Or should we have both options????
 *
 * @todo Use std::vector for building a list of pointers to observations
 * instead of a proper memory handling implementation.
 ***************************************************************************/
class GObservations {

    // I/O friends
    friend std::ostream& operator<<(std::ostream& os, const GObservations& obs);
    friend GLog&         operator<<(GLog& log,        const GObservations& obs);

public:
    // Constructors and destructors
    GObservations(void);
    GObservations(const GObservations& obs);
    explicit GObservations(const std::string& filename);
    virtual ~GObservations(void);

    // Operators
    GObservations&      operator= (const GObservations& obs);
    GObservation&       operator[](const int& index);
    const GObservation& operator[](const int& index) const;

    // Methods
    void        clear(void);
    int         size(void) const { return m_obs.size(); }
    void        append(GObservation& obs);
    void        load(const std::string& filename);
    void        save(const std::string& filename) const;
    void        read(const GXml& xml);
    void        write(GXml& xml) const;
    void        models(const GModels& models) { m_models=models;}
    void        models(const std::string& filename);
    GModels&    models(void) { return m_models; }
    void        optimize(GOptimizer& opt);
    double      npred(void) const { return m_npred; }
    std::string print(void) const;

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
        optimizer(GObservations* obs);
        optimizer(const optimizer& fct);
        ~optimizer(void);

        // Operators
        optimizer& operator= (const optimizer& fct);

        // Methods
        void           eval(const GOptimizerPars& pars);
        void           poisson_unbinned(const GObservation& obs, const GOptimizerPars& pars);
        void           poisson_unbinned(const GObservation& obs, const GOptimizerPars& pars, GSparseMatrix& covar, GVector& mgrad, double& value, GVector& gradient);
        void           poisson_binned(const GObservation& obs, const GOptimizerPars& pars);
        void           poisson_binned(const GObservation& obs, const GOptimizerPars& pars, GSparseMatrix& covar, GVector& mgrad, double& value, double& npred, GVector& gradient);
        void           gaussian_binned(const GObservation& obs, const GOptimizerPars& pars);
        void           gaussian_binned(const GObservation& obs, const GOptimizerPars& pars, GSparseMatrix& covar, GVector& mgrad, double& value, double& npred, GVector& gradient);
        double*        value(void) { return &m_value; }
        double         npred(void) const { return m_npred; }
        GVector*       gradient(void) { return m_gradient; }
        GSparseMatrix* covar(void) { return m_covar; }
    protected:
        void           init_members(void);
        void           copy_members(const optimizer& fct);
        void           free_members(void);
        double         m_value;       //!< Function value
        double         m_npred;       //!< Total number of predicted events
        double         m_minmod;      //!< Minimum model value
        double         m_minerr;      //!< Minimum error value
        GVector*       m_gradient;    //!< Pointer to gradient vector
        GSparseMatrix* m_covar;       //!< Pointer to covariance matrix
        GObservations* m_this;        //!< Pointer to GObservations object
        GVector*       m_wrk_grad;    //!< Pointer to working gradient vector
    };

protected:
    // Protected methods
    void           init_members(void);
    void           copy_members(const GObservations& obs);
    void           free_members(void);

    // Protected members
    std::vector<GObservation*> m_obs;      //!< List of observations
    GModels                    m_models;   //!< Models
    double                     m_npred;    //!< Total number of predicted events

};

#endif /* GOBSERVATIONS_HPP */
