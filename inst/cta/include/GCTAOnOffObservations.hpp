/***************************************************************************
 *    GCTAOnOffObservations.hpp - CTA on-off observation container class   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013 by Pierrick Martin                                  *
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
 * @file GCTAOnOffObservations.hpp
 * @brief CTA on-off observation container class definition
 * @author Pierrick Martin
 */

#ifndef GCTAONOFFOBSERVATIONS_HPP
#define GCTAONOFFOBSERVATIONS_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GContainer.hpp"
#include "GCTAOnOffObservation.hpp"
#include "GModels.hpp"


/***********************************************************************//**
 * @class GCTAOnOffObservations
 *
 * @brief ON/OFF Observation container class
 *
 * This class in a container of GCTAOnOffObservation objects.
 * Similar to a GObservations object, it holds a model container but for the
 * moment only the spectral part of it will be used. If no spectral model is
 * provided, the optimization cannot be performed
 *
 * The method optimize() performs the optimization of the spectral model
 * parameters applicable to all observations, together with the level of
 * background in each observation (which are nuisance parameters).
 * 
 ***************************************************************************/
class GCTAOnOffObservations : public GContainer {

public:
    // Constructors and destructors
    GCTAOnOffObservations(void);
    GCTAOnOffObservations(const GCTAOnOffObservations& obs);
    explicit GCTAOnOffObservations(const std::string& filename);
    virtual ~GCTAOnOffObservations(void);

    // Operators
    GCTAOnOffObservations&      operator=(const GCTAOnOffObservations& obs);
    GCTAOnOffObservation*       operator[](const int& index);
    const GCTAOnOffObservation* operator[](const int& index) const;

    // Methods
    void                        clear(void);
    GCTAOnOffObservations*      clone(void) const;
    GCTAOnOffObservation*       at(const int& index);
    const GCTAOnOffObservation* at(const int& index) const;
    int                         size(void) const;
    bool                        is_empty(void) const;
    GCTAOnOffObservation*       set(const int& index, const GCTAOnOffObservation& obs);
    GCTAOnOffObservation*       append(const GCTAOnOffObservation& obs);
    GCTAOnOffObservation*       insert(const int& index, const GCTAOnOffObservation& obs);
    void                        remove(const int& index);
    void                        reserve(const int& num);
    void                        extend(const GCTAOnOffObservations& obs);
    bool                        contains(const std::string& instrument,
                                         const std::string& id) const;
	void                        load(const std::string& filename);
    void                        save(const std::string& filename) const;
    void                        read(const GXml& xml);
    void                        write(GXml& xml) const;
	void                        models(const GModels& models);
    void                        models(const std::string& filename);
    const GModels&              models(void) const;	
    std::string                 print(const GChatter& chatter = NORMAL) const;
    void                        optimize(GOptimizer& opt);
	double                      npred(void) const;
	
	// Likelihood function
    class likelihood : public GOptimizerFunction {
    public:
        // Constructors and destructors
        likelihood(void);
        likelihood(GCTAOnOffObservation* obs);
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
        void set(GCTAOnOffObservations* obs);
        void eval(const GOptimizerPars& pars);
		
    protected:
        // Protected methods
        void           init_members(void);
        void           copy_members(const likelihood& fct);
        void           free_members(void);
		
        // Protected data members
        double                 m_value;       //!< Function value
        double                 m_npred;       //!< Total number of predicted events
        GVector*               m_gradient;    //!< Pointer to gradient vector
        GMatrixSparse*         m_curvature;   //!< Pointer to curvature matrix
        GCTAOnOffObservations* m_this;        //!< Pointer to GCTAOnOffObservations object
    };
	
    // Optimizer function access method
    const GCTAOnOffObservations::likelihood& function(void) const;
	
protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCTAOnOffObservations& obs);
    void free_members(void);
    int  get_index(const std::string& instrument,
                   const std::string& id) const;

    // Protected members
    std::vector<GCTAOnOffObservation*>  m_obs;    //!< List of observations
    GModels                             m_models; //!< List of models
	GCTAOnOffObservations::likelihood   m_fct;    //!< Optimizer function
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
GCTAOnOffObservation* GCTAOnOffObservations::operator[](const int& index)
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
const GCTAOnOffObservation* GCTAOnOffObservations::operator[](const int& index) const
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
int GCTAOnOffObservations::size(void) const
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
bool GCTAOnOffObservations::is_empty(void) const
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
void GCTAOnOffObservations::reserve(const int& num)
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
void GCTAOnOffObservations::models(const GModels& models)
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
const GModels& GCTAOnOffObservations::models(void) const
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
double GCTAOnOffObservations::npred(void) const
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
const GCTAOnOffObservations::likelihood& GCTAOnOffObservations::function(void) const
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
double GCTAOnOffObservations::likelihood::value(void)
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
double GCTAOnOffObservations::likelihood::npred(void) const
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
GVector* GCTAOnOffObservations::likelihood::gradient(void)
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
GMatrixSparse* GCTAOnOffObservations::likelihood::curvature(void)
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
void GCTAOnOffObservations::likelihood::set(GCTAOnOffObservations* obs)
{
    m_this = obs;
    return;
}


#endif /* GCTAONOFFOBSERVATIONS_HPP */
