/***************************************************************************
 *                     GModels.hpp - Model container class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2014 by Juergen Knoedlseder                         *
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
#include "GContainer.hpp"
#include "GModel.hpp"
#include "GOptimizerPars.hpp"
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
 * describe the data. The names of all models in the container have to be
 * unique, i.e. every name can occur only once. This allows for accessing the
 * models by name and by index.
 *
 * The GModels class provides methods to manage and to access the models
 * in the container. The number of models in the container is retrieved
 * using the size() method. The is_empty() method can be used to check
 * whether the container is empty or whether it contains models:
 *
 *     GModels models;                 // Allocate container
 *     int n = models.size();          // Number of models in container
 *     if (models.is_empty())           // Check for emptiness
 *         std::cout << "Empty container";
 *
 * Access operators exist for accessing of models by index or by name:
 *
 *     GModel* mptr = models[i];       // Get i'th model
 *     GModel* mptr = models["Crab"];  // Get a model named "Crab"
 *
 * The index access operator does not check the validity of the provided
 * index. For index validation, use the at() method:
 *
 *     GModel* mptr = models.at(i);    // Get i'th model with index check
 *
 * The append() method add a model to the container:
 *
 *     models.append(model);           // Append model
 *
 * The append() method clones the model that is passed as argument. The
 * method returns a pointer to the cloned model so that the attributes of
 * the cloned model can be manipulated:
 *
 *     GModel* mptr = models.append(model);
 *
 * The insert() methods insert a model before a given index or before
 * a model with a given name (the methods also return a pointer to the
 * cloned model):
 *
 *     models.insert(i, model);        // Insert before i'th model
 *     models.insert("Crab", model);   // Insert before "Crab" model
 *
 * The set() methods replace an existing model by index or by name (also
 * these methods return a pointer to the cloned model):
 *
 *     models.set(i, model);           // Replace i'th model
 *     models.set("Crab", model);      // Replace "Crab" model
 *
 * The remove() methods remove an existing model by index or by name:
 *
 *     models.remove(i);               // Remove i'th model
 *     models.remove("Crab");          // Remove "Crab" model
 *
 * The existence of a model with a given name can be checked using
 *
 *     if (models.contains("Crab"))
 *         std::cout << "We have the Crab!";
 *
 * The extend() method extends a container by all models found in another
 * container:
 *
 *     models.extend(other_models);    // Extend container
 *
 * For repeated container manipulations, a given @p number of model slots can
 * be reserved using
 *
 *     models.reserve(number);         // Reserves model slots
 *
 * which will speed up the memory allocations for new models.
 *
 * Models can be saved into or loaded from an XML file using
 *
 *     models.save("mymodels.xml");    // Save models in XML file
 *     models.load("mymodels.xml");    // Load models from XML file
 *
 * The models can also be loaded upon construction from an XML file:
 *
 *     GModels models("mymodels.xml"); // Construct by loading models from XML file
 *
 * The class can also directly operate on a GXml object using the read() and
 * write() methods:
 *
 *     GXml xml;                       // Allocate GXml object
 *     models.write(xml);              // Write into GXml object
 *     models.read(xml);               // Read models from GXml object
 *
 * The sum of all models in the container are evaluated for a given @p event
 * and @p observation using the eval() or eval_gradients() methods:
 *
 *     double value = models.eval(event, observation);
 *     double value = models.eval_gradients(event, observation);
 *
 * The eval_gradients() method sets the parameter gradients for all free
 * model parameters that have an analytical parameter gradient.
 *
 * The only member of GModels is a list of model pointers. The class handles
 * the proper allocation and deallocation of the model memory.
 ***************************************************************************/
class GModels : public GContainer {

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
    void           clear(void);
    GModels*       clone(void) const;
    std::string    classname(void) const;
    GModel*        at(const int& index);
    const GModel*  at(const int& index) const;
    int            size(void) const;
    bool           is_empty(void) const;
    GModel*        set(const int& index, const GModel& model);
    GModel*        set(const std::string& name, const GModel& model);
    GModel*        append(const GModel& model);
    GModel*        insert(const int& index, const GModel& model);
    GModel*        insert(const std::string& name, const GModel& model);
    void           remove(const int& index);
    void           remove(const std::string& name);
    void           reserve(const int& num);
    void           extend(const GModels& models);
    bool           contains(const std::string& name) const;
    void           load(const std::string& filename);
    void           save(const std::string& filename) const;
    void           read(const GXml& xml);
    void           write(GXml& xml) const;
    int            npars(void) const;
    GOptimizerPars pars(void);
    double         eval(const GEvent& event, const GObservation& obs) const;
    double         eval_gradients(const GEvent& event, const GObservation& obs) const;
    std::string    print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void          init_members(void);
    void          copy_members(const GModels& models);
    void          free_members(void);
    int           get_index(const std::string& name) const;

    // Proteced members
    std::vector<GModel*> m_models;  //!< List of models
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GModels").
 ***************************************************************************/
inline
std::string GModels::classname(void) const
{
    return ("GModels");
}


/***********************************************************************//**
 * @brief Return pointer to model
 *
 * @param[in] index Model index [0,...,size()-1].
 *
 * Returns a pointer to the model with the specified @p index.
 ***************************************************************************/
inline
GModel* GModels::operator[](const int& index)
{
    return (m_models[index]);
}


/***********************************************************************//**
 * @brief Return pointer to model (const version)
 *
 * @param[in] index Model index [0,...,size()-1].
 *
 * Returns a const pointer to the model with the specified @p index.
 ***************************************************************************/
inline
const GModel* GModels::operator[](const int& index) const
{
    return (m_models[index]);
}


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
bool GModels::is_empty(void) const
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
