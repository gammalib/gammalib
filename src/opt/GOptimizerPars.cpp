/***************************************************************************
 *         GOptimizerPars.cpp - Optimizer parameter container class        *
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
 * @file GOptimizerPars.cpp
 * @brief Optimizer parameter container class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GOptimizerPars.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_ACCESS                   "GOptimizerPars::operator[](std::string&)"
#define G_AT1                                      "GOptimizerPars::at(int&)"
#define G_AT2                              "GOptimizerPars::at(std::string&)"
#define G_SET1                    "GOptimizerPars::set(int&, GOptimizerPar&)"
#define G_SET2            "GOptimizerPars::set(std::string&, GOptimizerPar&)"
#define G_ATTACH1              "GOptimizerPars::attach(int&, GOptimizerPar*)"
#define G_ATTACH2      "GOptimizerPars::attach(std::string&, GOptimizerPar*)"
#define G_INSERT1              "GOptimizerPars::insert(int&, GOptimizerPar&)"
#define G_INSERT2      "GOptimizerPars::insert(std::string&, GOptimizerPar&)"
#define G_REMOVE1                              "GOptimizerPars::remove(int&)"
#define G_REMOVE2                      "GOptimizerPars::remove(std::string&)"


/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GOptimizerPars::GOptimizerPars(void)
{
    // Initialise members
    init_members();
  
    // Return
    return;
}


/***********************************************************************//**
 * @brief Parameter constructor
 *
 * @param[in] number Number of parameters to allocate.
 *
 * Allocates @p number parameters in the parameter container.
 ***************************************************************************/
GOptimizerPars::GOptimizerPars(const int& number)
{
    // Initialise members
    init_members();

    // Allocate parameters
    for (int i = 0; i < number; ++i) {
        m_alloc.push_back(true);
        m_pars.push_back(new GOptimizerPar);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] pars Optimizer parameters.
 ***************************************************************************/
GOptimizerPars::GOptimizerPars(const GOptimizerPars& pars)
{ 
    // Initialise members
    init_members();

    // Copy members
    copy_members(pars);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GOptimizerPars::~GOptimizerPars()
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Operators                                =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] pars Optimizer parameters.
 * @return Optimizer parameters.
 ***************************************************************************/
GOptimizerPars& GOptimizerPars::operator=(const GOptimizerPars& pars)
{ 
    // Execute only if object is not identical
    if (this != &pars) {

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(pars);

    } // endif: object was not identical
  
    // Return
    return *this;
}


/***********************************************************************//**
 * @brief Return pointer to parameter
 *
 * @param[in] name Parameter name.
 *
 * @exception GException::invalid_argument
 *            Parameter with specified @p name not found in container.
 *
 * Returns a pointer to the parameter with the specified @p name.
 ***************************************************************************/
GOptimizerPar* GOptimizerPars::operator[](const std::string& name)
{
    // Get parameter index
    int index = get_index(name);

    // Throw exception if model name was not found
    if (index == -1) {
        std::string msg = "Parameter \""+name+"\" not found in parameter"
                          " container.\nPlease specify a valid parameter"
                          " name.";
        throw GException::invalid_argument(G_ACCESS, msg);
    }

    // Return pointer
    return m_pars[index];
}


/***********************************************************************//**
 * @brief Return pointer to parameter (const version)
 *
 * @param[in] name Parameter name.
 *
 * @exception GException::invalid_argument
 *            Parameter with specified @p name not found in container.
 *
 * Returns a pointer to the parameter with the specified @p name.
 ***************************************************************************/
const GOptimizerPar* GOptimizerPars::operator[](const std::string& name) const
{
    // Get parameter index
    int index = get_index(name);

    // Throw exception if model name was not found
    if (index == -1) {
        std::string msg = "Parameter \""+name+"\" not found in parameter"
                          " container.\nPlease specify a valid parameter"
                          " name.";
        throw GException::invalid_argument(G_ACCESS, msg);
    }

    // Return pointer
    return m_pars[index];
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear parameter container
 *
 * Removes all parameters from the container.
 ***************************************************************************/
void GOptimizerPars::clear(void)
{
    // Free class members
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone parameter container
 *
 * @return Pointer to deep copy of parameter container
 ***************************************************************************/
GOptimizerPars* GOptimizerPars::clone(void) const
{
    return new GOptimizerPars(*this);
}


/***********************************************************************//**
 * @brief Return pointer to parameter
 *
 * @param[in] index Parameter index [0,...,size()-1].
 * @return Pointer to parameter.
 *
 * @exception GException::out_of_range
 *            Extension number is out of range.
 *
 * Returns a pointer to the parameter with the specified @p index.
 ***************************************************************************/
GOptimizerPar* GOptimizerPars::at(const int& index)
{
    // Compile option: raise an exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT1, "Parameter index", index, size());
    }
    #endif

    // Return pointer
    return (m_pars[index]);
}


/***********************************************************************//**
 * @brief Return pointer to parameter (const version)
 *
 * @param[in] index Parameter index [0,...,size()-1].
 * @return Pointer to parameter.
 *
 * @exception GException::out_of_range
 *            Extension number is out of range.
 *
 * Returns a pointer to the parameter with the specified @p index.
 ***************************************************************************/
const GOptimizerPar* GOptimizerPars::at(const int& index) const
{
    // Compile option: raise an exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT1, "Parameter index", index, size());
    }
    #endif

    // Return pointer
    return (m_pars[index]);
}


/***********************************************************************//**
 * @brief Return number of free parameters
 *
 * @return Number of free parameters.
 *
 * Determines the number of free parameters by collecting statistics from all
 * model parameters.
 ***************************************************************************/
int GOptimizerPars::nfree(void) const
{
    // Initialise number of free parameters
    int nfree = 0;
    
    // Collect all free parameters
    for (int i = 0; i < size(); ++i) {
        if (m_pars[i]->isfree()) {
            nfree++;
        }
    }
    
    // Return
    return nfree;
}


/***********************************************************************//**
 * @brief Set parameter in container
 *
 * @param[in] index Parameter index [0,...,size()-1].
 * @param[in] par Parameter.
 * @return Pointer to deep copy of parameter.
 *
 * @exception GException::out_of_range
 *            Parameter index is out of range.
 *
 * Set parameter in the container. A deep copy of the parameter will be made.
 ***************************************************************************/
GOptimizerPar* GOptimizerPars::set(const int& index, const GOptimizerPar& par)
{
    // Compile option: raise an exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_SET1, "Parameter index", index, size());
    }
    #endif

    // Delete parameter if it has been allocated
    if (m_alloc[index]) {
        delete m_pars[index];
    }

    // Clone parameter
    GOptimizerPar* ptr = par.clone();

    // Set parameter
    m_alloc[index] = true;
    m_pars[index]  = ptr;

    // Return parameter pointer
    return ptr;
}


/***********************************************************************//**
 * @brief Set parameter in container
 *
 * @param[in] name Parameter name.
 * @param[in] par Parameter.
 * @return Pointer to deep copy of parameter.
 *
 * @exception GException::invalid_argument
 *            Specified parameter @p name not found in container.
 ***************************************************************************/
GOptimizerPar* GOptimizerPars::set(const std::string& name, const GOptimizerPar& par)
{
    // Get parameter index
    int index = get_index(name);

    // Throw exception if model name was not found
    if (index == -1) {
        std::string msg = "Parameter \""+name+"\" not found in parameter"
                          " container.\nPlease specify a valid parameter"
                          " name.";
        throw GException::invalid_argument(G_SET2, msg);
    }

    // Set parameter
    return set(index, par);
}


/***********************************************************************//**
 * @brief Attach parameter to container
 *
 * @param[in] par Parameter pointer.
 ***************************************************************************/
void GOptimizerPars::attach(GOptimizerPar* par)
{
    // Push pointer in vector if it is valid
    if (par != NULL) {
        m_alloc.push_back(false);
        m_pars.push_back(par);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Attach parameter to container at the specified index
 *
 * @param[in] index Parameter index [0,...,size()-1].
 * @param[in] par Parameter pointer.
 *
 * @exception GException::out_of_range
 *            Parameter @p index is not valid.
 *
 * Attaches parameter pointer at the specified @p index.
 ***************************************************************************/
void GOptimizerPars::attach(const int& index, GOptimizerPar* par)
{
    // Compile option: raise an exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_ATTACH1, "Parameter index", index, size());
    }
    #endif

    // Delete parameter if it has been allocated
    if (m_alloc[index]) {
        delete m_pars[index];
    }

    // Set parameter
    m_alloc[index] = false;
    m_pars[index]  = par;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Attach parameter to container at the specified index
 *
 * @param[in] name Parameter name.
 * @param[in] par Parameter pointer.
 *
 * @exception GException::out_of_range
 *            Parameter @p index is not valid.
 *
 * Attaches parameter pointer at the specified @p index.
 ***************************************************************************/
void GOptimizerPars::attach(const std::string& name, GOptimizerPar* par)
{
    // Get parameter index
    int index = get_index(name);

    // Throw exception if model name was not found
    if (index == -1) {
        std::string msg = "Parameter \""+name+"\" not found in parameter"
                          " container.\nPlease specify a valid parameter"
                          " name.";
        throw GException::invalid_argument(G_ATTACH2, msg);
    }

    // Remove parameter
    attach(index, par);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Insert parameter into container
 *
 * @param[in] index Parameter index [0,...,size()-1].
 * @param[in] par Parameter.
 * @return Pointer to deep copy of parameter.
 *
 * @exception GException::out_of_range
 *            Parameter index is out of range.
 *
 * Set parameter in the container. A deep copy of the parameter will be made.
 ***************************************************************************/
GOptimizerPar* GOptimizerPars::insert(const int& index, const GOptimizerPar& par)
{
    // Compile option: raise an exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_INSERT1, "Parameter index", index, size());
    }
    #endif

    // Clone parameter
    GOptimizerPar* ptr = par.clone();

    // Inserts deep copy of parameter
    m_alloc.insert(m_alloc.begin()+index, true);
    m_pars.insert(m_pars.begin()+index, ptr);

    // Return parameter pointer
    return ptr;
}


/***********************************************************************//**
 * @brief Insert parameter into container
 *
 * @param[in] name Parameter name.
 * @param[in] par Parameter.
 * @return Pointer to deep copy of parameter.
 *
 * @exception GException::invalid_argument
 *            Specified parameter @p name not found in container.
 ***************************************************************************/
GOptimizerPar* GOptimizerPars::insert(const std::string& name,
                                      const GOptimizerPar& par)
{
    // Get parameter index
    int index = get_index(name);

    // Throw exception if model name was not found
    if (index == -1) {
        std::string msg = "Parameter \""+name+"\" not found in parameter"
                          " container.\nPlease specify a valid parameter"
                          " name.";
        throw GException::invalid_argument(G_INSERT2, msg);
    }

    // Insert parameter
    return insert(index, par);
}


/***********************************************************************//**
 * @brief Remove parameter from container
 *
 * @param[in] index Parameter index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Parameter index is out of range.
 *
 * Remove parameter of specified @p index from container.
 ***************************************************************************/
void GOptimizerPars::remove(const int& index)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_REMOVE1, "Parameter index",
                                       index, size());
    }
    #endif

    // Delete parameter if it has been allocated
    if (m_alloc[index]) {
        delete m_pars[index];
    }

    // Erase parameter from container
    m_alloc.erase(m_alloc.begin() + index);
    m_pars.erase(m_pars.begin() + index);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Remove parameter from container
 *
 * @param[in] name Parameter name.
 *
 * @exception GException::invalid_argument
 *            Specified parameter @p name not found in container.
 ***************************************************************************/
void GOptimizerPars::remove(const std::string& name)
{
    // Get parameter index
    int index = get_index(name);

    // Throw exception if model name was not found
    if (index == -1) {
        std::string msg = "Parameter \""+name+"\" not found in parameter"
                          " container.\nPlease specify a valid parameter"
                          " name.";
        throw GException::invalid_argument(G_REMOVE2, msg);
    }

    // Remove parameter
    remove(index);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Append parameter container
 *
 * @param[in] pars Parameter container.
 *
 * Append parameter container to the container. All parameters that were
 * allocated within @p pars will be clones, for the other parameters, the
 * pointers will be copied.
 ***************************************************************************/
void GOptimizerPars::extend(const GOptimizerPars& pars)
{
    // Do nothing if parameter container is empty
    if (!pars.isempty()) {

        // Get size. Note that we extract the size first to avoid an
        // endless loop that arises when a container is appended to
        // itself.
        int num = pars.size();

        // Reserve enough space
        reserve(size() + num);

        // Loop over all parameters 
        for (int i = 0; i < num; ++i) {
            GOptimizerPar* par = (pars.m_alloc[i]) ? pars.m_pars[i]->clone()
                                                   : pars.m_pars[i];
            m_pars.push_back(par);
        }

    } // endif: parameter container was not empty
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Signals if parameter name exists
 *
 * @param[in] name Parameter name.
 * @return True if parameter with specified @p name exists.
 *
 * Searches all parameter names for a match with the specified @p name. If
 * the specified name has been found, true is returned.
 ***************************************************************************/
bool GOptimizerPars::contains(const std::string& name) const
{
    // Get model index
    int index = get_index(name);

    // Return
    return (index != -1);
}


/***********************************************************************//**
 * @brief Print parameters
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing parameter container information.
 *
 * Prints all parameters into a string.
 ***************************************************************************/
std::string GOptimizerPars::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GOptimizerPars ===");

        // Append information
        result.append("\n"+gammalib::parformat("Number of parameters"));
        result.append(gammalib::str(size()));

        // Append parameters
        for (int i = 0; i < size(); ++i) {
            result.append("\n"+m_pars[i]->print(chatter));
        }

    } // endif: chatter was not silent

    // Return result
    return result;
}



/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GOptimizerPars::init_members(void)
{
    // Initialise members
    m_alloc.clear();
    m_pars.clear();
  
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] pars Optimizer parameters.
 ***************************************************************************/
void GOptimizerPars::copy_members(const GOptimizerPars& pars)
{
    // Copy members
    m_alloc = pars.m_alloc;

    // Clone or copy parameter pointers, depending on whether they have
    // been allocated or not in the instance from which we copy
    for (int i = 0; i < size(); ++i) {
        GOptimizerPar* par = (pars.m_alloc[i]) ? pars.m_pars[i]->clone()
                                               : pars.m_pars[i];
        m_pars.push_back(par);
    }
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GOptimizerPars::free_members(void)
{
    // Free all allocated parameters
    for (int i = 0; i < size(); ++i) {
        if (m_alloc[i] && m_pars[i] != NULL) {
            delete m_pars[i];
            m_pars[i] = NULL;
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return parameter index by name
 *
 * @param[in] name Parameter name.
 * @return Parameter index (-1 if parameter name has not been found)
 *
 * Returns parameter index based on the specified @p name. If no parameter
 * with the specified @p name is found the method returns -1.
 ***************************************************************************/
int GOptimizerPars::get_index(const std::string& name) const
{
    // Initialise index
    int index = -1;

    // Search parameter with specified name
    for (int i = 0; i < size(); ++i) {
        if (m_pars[i]->name() == name) {
            index = i;
            break;
        }
    }

    // Return index
    return index;
}
