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
 * gamma-ray observations and manages their scientific analysis. For a given
 * instruments, the identifiers of all observations in the container have to
 * be unique, i.e. every identifier can occur only once. This allows for
 * accessing the observations by instrument and by identifier.
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
 * returning a reference to the observation. The operator does not perform
 * index range checking. If range checking is required, use the at() method.
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
 * npred() method returns the total number of events that are predicted
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
    GObservations&      operator=(const GObservations& obs);
    GObservation*       operator[](const int& index);
    const GObservation* operator[](const int& index) const;

    // Methods
    void                clear(void);
    GObservations*      clone(void) const;
    GObservation*       at(const int& index);
    const GObservation* at(const int& index) const;
    int                 size(void) const;
    bool                is_empty(void) const;
    GObservation*       set(const int& index, const GObservation& obs);
    GObservation*       append(const GObservation& obs);
    GObservation*       insert(const int& index, const GObservation& obs);
    void                remove(const int& index);
    void                reserve(const int& num);
    void                extend(const GObservations& obs);
    bool                contains(const std::string& instrument,
                                 const std::string& id) const;
    void                load(const std::string& filename);
    void                save(const std::string& filename) const;
    void                read(const GXml& xml);
    void                write(GXml& xml) const;
    void                models(const GModels& models);
    void                models(const std::string& filename);
    const GModels&      models(void) const;
    void                optimize(GOptimizer& opt);
    double              npred(void) const;
    std::string         print(const GChatter& chatter = NORMAL) const;

    // Likelihood function
    class likelihood : public GOptimizerFunction {
    public:
        // Constructors and destructors
        likelihood(void);
        likelihood(GObservations* obs);
        likelihood(const likelihood& fct);
        ~likelihood(void);

        // Operators
        likelihood& operator=(const likelihood& fct);

        // Implemented pure virtual base class methods
        double         value(void);
        double         npred(void) const;
        GVector*       gradient(void);
        GMatrixSparse* curvature(void);

        // Other methods
        void set(GObservations* obs);
        void eval(const GOptimizerPars& pars);

    protected:
        // Protected methods
        void           init_members(void);
        void           copy_members(const likelihood& fct);
        void           free_members(void);

        // Protected data members
        double         m_value;       //!< Function value
        double         m_npred;       //!< Total number of predicted events
        GVector*       m_gradient;    //!< Pointer to gradient vector
        GMatrixSparse* m_curvature;   //!< Pointer to curvature matrix
        GObservations* m_this;        //!< Pointer to GObservations object
    };

    // Optimizer function access method
    const GObservations::likelihood& function(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GObservations& obs);
    void free_members(void);
    int  get_index(const std::string& instrument,
                   const std::string& id) const;

    // Protected members
    std::vector<GObservation*> m_obs;    //!< List of observations
    GModels                    m_models; //!< List of models
    GObservations::likelihood  m_fct;    //!< Optimizer function
};


/***********************************************************************//**
 * @brief Return pointer to observation
 *
 * @param[in] index Observation index [0,...,size()-1].
 * @return Observation.
 *
 * Returns a pointer to the observation with specified @p index.
 ***************************************************************************/
inline
GObservation* GObservations::operator[](const int& index)
{
    return m_obs[index];
}


/***********************************************************************//**
 * @brief Return pointer to observation (const version)
 *
 * @param[in] index Observation index [0,...,size()-1].
 *
 * Returns a const pointer to the observation with specified @p index.
 ***************************************************************************/
inline
const GObservation* GObservations::operator[](const int& index) const
{
    return m_obs[index];
}


/***********************************************************************//**
 * @brief Return number of observations in container
 *
 * @return Number of observations in container.
 *
 * Returns the number of observations in the observation container.
 ***************************************************************************/
inline
int GObservations::size(void) const
{
    return (m_obs.size());
}


/***********************************************************************//**
 * @brief Signals if there are no observations in container
 *
 * @return True if container is empty, false otherwise.
 *
 * Signals if the observation container does not contain any observation.
 ***************************************************************************/
inline
bool GObservations::is_empty(void) const
{
    return (m_obs.empty());
}


/***********************************************************************//**
 * @brief Reserves space for observations in container
 *
 * @param[in] num Number of observations.
 *
 * Reserves space for @p num observations in the container.
 ***************************************************************************/
inline
void GObservations::reserve(const int& num)
{
    m_obs.reserve(num);
    return;
}


/***********************************************************************//**
 * @brief Set model container
 *
 * @param[in] models Model container.
 *
 * Sets the model container for the observations.
 ***************************************************************************/
inline
void GObservations::models(const GModels& models)
{
    m_models = models;
    return;
}


/***********************************************************************//**
 * @brief Return model container
 *
 * @return Model container.
 *
 * Returns the model container of the observations.
 ***************************************************************************/
inline
const GModels& GObservations::models(void) const
{
    return m_models;
}


/***********************************************************************//**
 * @brief Return total number of predicted events after model fitting
 *
 * @return Total number of predicted events after model fitting.
 *
 * Returns the total number of events that is predicted by the models after
 * they have been fitted to the data. This number if computed following a
 * call of the optimize() method. If this method has never been called, the
 * npred() method returns 0.
 ***************************************************************************/
inline
double GObservations::npred(void) const
{
    return (m_fct.npred());
}


/***********************************************************************//**
 * @brief Return likelihood function
 *
 * @return Reference to likelihood function.
 *
 * Returns a reference to the likelihood function.
 ***************************************************************************/
inline
const GObservations::likelihood& GObservations::function(void) const
{
    return m_fct;
}


/***********************************************************************//**
 * @brief Return likelihood function value
 *
 * @return Likelihood function value.
 *
 * Returns the actual function value of the likelihood function.
 ***************************************************************************/
inline
double GObservations::likelihood::value(void)
{
    return m_value;
}


/***********************************************************************//**
 * @brief Return total number of predicted events
 *
 * @return Total number of predicted events.
 *
 * Returns the total number of events that is predicted by the models after
 * they have been fitted to the data.
 ***************************************************************************/
inline
double GObservations::likelihood::npred(void) const
{
    return m_npred;
}


/***********************************************************************//**
 * @brief Return pointer to gradient vector
 *
 * @return Pointer to gradient vector.
 *
 * Returns a pointer to the parameter gradient vector.
 ***************************************************************************/
inline
GVector* GObservations::likelihood::gradient(void)
{
    return m_gradient;
}


/***********************************************************************//**
 * @brief Return pointer to curvature matrix
 *
 * @return Pointer to curvature matrix.
 *
 * Returns a pointer to the parameter curvature matrix.
 ***************************************************************************/
inline
GMatrixSparse* GObservations::likelihood::curvature(void)
{
    return m_curvature;
}


/***********************************************************************//**
 * @brief Set observation container
 *
 * @param[in] obs Pointer to observation container.
 *
 * Sets the pointer to the observation container for which the optimizer
 * class should be used.
 ***************************************************************************/
inline
void GObservations::likelihood::set(GObservations* obs)
{
    m_this = obs;
    return;
}

#endif /* GOBSERVATIONS_HPP */
