/***************************************************************************
 *                     GModels.hpp - Model container class                 *
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
 * @file GModels.hpp
 * @brief Model container class definition
 * @author Juergen Knoedlseder
 */

#ifndef GMODELS_HPP
#define GMODELS_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GOptimizerPars.hpp"
#include "GModel.hpp"
#include "GXml.hpp"

/* __ Forward declarations _______________________________________________ */
class GEvent;
class GObservation;


/***********************************************************************//**
 * @class GModels
 *
 * @brief Model container class
 *
 * This container class collects models of type GModel that are used to
 * describe the gamma-ray data. Each model has a number of parameters
 * that are implemented using the GModelPar class. This container class
 * provides methods to manage the container class, to access and to evaluate
 * the models.
 *
 * The only member of this class is a list of model pointers. GModels handles
 * the proper allocation and deallocation of the model memory.
 *
 * GModels derives from GOptimizerPars which contains a flat array of
 * model parameters. This flat array is set using the protected
 * set_pointers() method. This method is called after each manipulation of
 * the list of model pointers to ensure consistency between the models and
 * the flat array of parameter pointers. The optimizer function
 * GOptimizerFunction will manipulate the parameters in that flat array,
 * and as this array points towards the parameters stored in GModels, it will
 * directly modify the model parameters of the container class.
 ***************************************************************************/
class GModels : public GOptimizerPars {

public:
    // Constructors and destructors
    GModels(void);
    GModels(const GModels& models);
    explicit GModels(const std::string& filename);
    virtual ~GModels(void);

    // Operators
    GModels&      operator=(const GModels& models);
    GModel*       operator[](const int& index);
    const GModel* operator[](const int& index) const;
    GModel*       operator[](const std::string& name);
    const GModel* operator[](const std::string& name) const;

    // Methods
    void          clear(void);
    GModels*      clone(void) const;
    GModel*       at(const int& index);
    const GModel* at(const int& index) const;
    int           size(void) const;
    bool          isempty(void) const;
    GModel*       set(const int& index, const GModel& model);
    GModel*       set(const std::string& name, const GModel& model);
    GModel*       append(const GModel& model);
    GModel*       insert(const int& index, const GModel& model);
    GModel*       insert(const std::string& name, const GModel& model);
    void          remove(const int& index);
    void          remove(const std::string& name);
    void          reserve(const int& num);
    void          extend(const GModels& models);
    bool          hasmodel(const std::string& name) const;
    void          load(const std::string& filename);
    void          save(const std::string& filename) const;
    void          read(const GXml& xml);
    void          write(GXml& xml) const;
    double        eval(const GEvent& event, const GObservation& obs) const;
    double        eval_gradients(const GEvent& event, const GObservation& obs) const;
    std::string   print(void) const;

protected:
    // Protected methods
    void          init_members(void);
    void          copy_members(const GModels& models);
    void          free_members(void);
    void          set_pointers(void);
    int           get_index(const std::string& name) const;

    // Proteced members
    std::vector<GModel*> m_models;  //!< List of models
};


/***********************************************************************//**
 * @brief Return number of models in container
 *
 * @return Number of models in container.
 *
 * Returns the number of models in the model container.
 ***************************************************************************/
inline
int GModels::size(void) const
{
    return (m_models.size());
}


/***********************************************************************//**
 * @brief Signals if there are no models in container
 *
 * @return True if container is empty, false otherwise.
 *
 * Signals if the model container does not contain any model.
 ***************************************************************************/
inline
bool GModels::isempty(void) const
{
    return (m_models.empty());
}


/***********************************************************************//**
 * @brief Reserves space for models in container
 *
 * @param[in] num Number of models
 *
 * Reserves space for @p num models in the container.
 ***************************************************************************/
inline
void GModels::reserve(const int& num)
{
    m_models.reserve(num);
    return;
}

#endif /* GMODELS_HPP */
