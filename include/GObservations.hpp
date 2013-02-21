/***************************************************************************
 *              GObservations.hpp - Observation container class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2013 by Juergen Knoedlseder                         *
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
 * @author Juergen Knoedlseder
 */

#ifndef GOBSERVATIONS_HPP
#define GOBSERVATIONS_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GContainer.hpp"
#include "GObservation.hpp"
#include "GEvent.hpp"
#include "GOptimizer.hpp"
#include "GOptimizerFunction.hpp"
#include "GModels.hpp"


/***********************************************************************//**
 * @class GObservations
 *
 * @brief Observation container class
 *
 * This is the main user interface class that provides high-level access to
 * gamma-ray observations and manages their scientific analysis.
 *
 * GObservations holds a list of pointers to GObservation objects which
 * implement gamma-ray observations. The class derives from the abstract
 * container interface class GContainer, and provides methods to manage the
 * list of observations.
 * The append() appends an observation to the container, the insert() method
 * inserts an observation at a specific position in the list (though the
 * order of the observations in the list is not relevant for the scientific
 * analysis).
 * The remove() method removes a specific observation from the list and the
 * extend() method extends the list of observations by appending another
 * list of observations to it.
 * The size() method provides the number of observations in the list.
 *
 * The operator[] provides access to an observation in the list by index,
 * returning a reference to the observation.
 *
 * The observation information may be stored into a XML file. The save()
 * method saves all information into a XML file while the load() method
 * loads the information into the list. Alternatively, a XML file
 * constructor can be used to construct an instance of GObservations by
 * loading observation information from a XML file. In addition, the
 * read() and write() method handle the reading and writing of observation
 * information from and into a XML object of type GXml.
 *
 * The class also holds a list of models that is implemented by the GModels
 * class. The models() method allows setting and retrieving the list of
 * models, as well as loading of models from an XML file.
 *
 * Based on the list of models, the optimize() method will optimize all
 * model parameters that are marked as free in the list of models. The
 * npred() method returns the total number of events that are prediced
 * by all models after optimization.
 *
 * GObservations also provides an optimizer class that is derived from
 * the abstract GOptimizerFunction base class. The GObservations::optimizer
 * class is the object that is used for model parameter optimization.
 ***************************************************************************/
class GObservations : public GContainer {

public:
    // Constructors and destructors
    GObservations(void);
    GObservations(const GObservations& obs);
    explicit GObservations(const std::string& filename);
    virtual ~GObservations(void);

    // Operators
    GObservations&      operator= (const GObservations& obs);
    GObservation*       operator[](const int& index);
    const GObservation* operator[](const int& index) const;

    // Methods
    void           clear(void);
    GObservations* clone(void) const;
    int            size(void) const { return m_obs.size(); }
    bool           isempty(void) const { return m_obs.empty(); }
    void           set(const int& index, const GObservation& obs);
    void           append(const GObservation& obs);
    void           insert(const int& index, const GObservation& obs);
    void           remove(const int& index);
    void           reserve(const int& num) { return m_obs.reserve(num); }
    void           extend(const GObservations& obs);
    void           load(const std::string& filename);
    void           save(const std::string& filename) const;
    void           read(const GXml& xml);
    void           write(GXml& xml) const;
    void           models(const GModels& models) { m_models=models;} //!< @brief Set model container
    void           models(const std::string& filename);
    const GModels& models(void) { return m_models; } //!< @brief Return model container
    void           optimize(GOptimizer& opt);
    double         npred(void) const { return m_fct.npred(); } //!< @brief Return total number of predicted events
    std::string    print(void) const;


    // Optimizer
    class optimizer : public GOptimizerFunction {
    public:
        // Constructors and destructors
        optimizer(void);
        optimizer(GObservations* obs);
        optimizer(const optimizer& fct);
        ~optimizer(void);

        // Operators
        optimizer& operator=(const optimizer& fct);

        // Implemented pure virtual base class methods
        double         value(void) { return m_value; }       //!< @brief Return optimizer function value
        double         npred(void) const { return m_npred; } //!< @brief Return predicted number of events
        GVector*       gradient(void) { return m_gradient; } //!< @brief Return pointer to gradient vector
        GSparseMatrix* covar(void) { return m_covar; }       //!< @brief Return pointer to covariance matrix

        // Other methods
        void set(GObservations* obs) { m_this=obs; } //!< @brief Set GObservations pointer
        void eval(const GOptimizerPars& pars);
        void poisson_unbinned(const GObservation&   obs,
                              const GOptimizerPars& pars);
        void poisson_unbinned(const GObservation&   obs,
                              const GOptimizerPars& pars,
                              GSparseMatrix&        covar,
                              GVector&              mgrad,
                              double&               value,
                              GVector&              gradient);
        void poisson_binned(const GObservation&   obs,
                            const GOptimizerPars& pars);
        void poisson_binned(const GObservation&   obs,
                            const GOptimizerPars& pars,
                            GSparseMatrix&        covar,
                            GVector&              mgrad,
                            double&               value,
                            double&               npred,
                            GVector&              gradient);
        void gaussian_binned(const GObservation&   obs,
                             const GOptimizerPars& pars);
        void gaussian_binned(const GObservation&   obs,
                             const GOptimizerPars& pars,
                             GSparseMatrix&        covar,
                             GVector&              mgrad,
                             double&               value,
                             double&               npred,
                             GVector&              gradient);

    protected:
        // Protected methods
        void           init_members(void);
        void           copy_members(const optimizer& fct);
        void           free_members(void);

        // Protected data members
        double         m_value;       //!< Function value
        double         m_npred;       //!< Total number of predicted events
        double         m_minmod;      //!< Minimum model value
        double         m_minerr;      //!< Minimum error value
        GVector*       m_gradient;    //!< Pointer to gradient vector
        GSparseMatrix* m_covar;       //!< Pointer to covariance matrix
        GVector*       m_wrk_grad;    //!< Pointer to working gradient vector
        GObservations* m_this;        //!< Pointer to GObservations object
    };

    // Optimizer access method
    const GObservations::optimizer& function(void) const { return m_fct; } //!< @brief Return optimizer function

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GObservations& obs);
    void free_members(void);

    // Protected members
    std::vector<GObservation*> m_obs;    //!< List of observations
    GModels                    m_models; //!< List of models
    GObservations::optimizer   m_fct;    //!< Optimizer function
};

#endif /* GOBSERVATIONS_HPP */
