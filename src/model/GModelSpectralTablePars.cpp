/***************************************************************************
 * GModelSpectralTablePars.cpp - Spectral table model par container class  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2019 by Juergen Knoedlseder                              *
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
 * @file GModelSpectralTablePars.cpp
 * @brief Spectral table model parameter container class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include "GException.hpp"
#include "GTools.hpp"
#include "GModelSpectralTablePars.hpp"

/* __ Constants __________________________________________________________ */

/* __ Globals ____________________________________________________________ */

/* __ Method name definitions ____________________________________________ */
#define G_ACCESS          "GModelSpectralTablePars::operator[](std::string&)"
#define G_AT                              "GModelSpectralTablePars::at(int&)"
#define G_SET1  "GModelSpectralTablePars::set(int&, GModelSpectralTablePar&)"
#define G_SET2                  "GModelSpectralTablePars::set(std::string&, "\
                                                   "GModelSpectralTablePar&)"
#define G_APPEND   "GModelSpectralTablePars::append(GModelSpectralTablePar&)"
#define G_INSERT1                    "GModelSpectralTablePars::insert(int&, "\
                                                   "GModelSpectralTablePar&)"
#define G_INSERT2            "GModelSpectralTablePars::insert(std::string&, "\
                                                   "GModelSpectralTablePar&)"
#define G_REMOVE1                     "GModelSpectralTablePars::remove(int&)"
#define G_REMOVE2  "GModelSpectralTablePars::remove(const std::string& name)"
#define G_EXTEND  "GModelSpectralTablePars::extend(GModelSpectralTablePars&)"

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
GModelSpectralTablePars::GModelSpectralTablePars(void)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] pars Table model parameters.
 ***************************************************************************/
GModelSpectralTablePars::GModelSpectralTablePars(const GModelSpectralTablePars& pars)
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
GModelSpectralTablePars::~GModelSpectralTablePars(void)
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
 * @param[in] pars Table model parameters.
 * @return Table model parameters.
 ***************************************************************************/
GModelSpectralTablePars& GModelSpectralTablePars::operator=(const GModelSpectralTablePars& pars)
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
 * @brief Return pointer to table model parameter (const version)
 *
 * @param[in] name Table model parameter name.
 *
 * @exception GException::invalid_argument
 *            Table model parameter with specified name not found in
 *            container.
 *
 * Returns a const pointer to the table model parameter with the specified
 * @p name.
 ***************************************************************************/
const GModelSpectralTablePar* GModelSpectralTablePars::operator[](const std::string& name) const
{
    // Get model index
    int index = get_index(name);

    // Throw exception if model name was not found
    if (index == -1) {
        std::string msg = "Table model parameter \""+name+"\" not found. "
                          "Please specify a valid parameter name.";
        throw GException::invalid_argument(G_ACCESS, msg);
    }

    // Return pointer
    return m_pars[index];
}


/*==========================================================================
 =                                                                         =
 =                              Public methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear table model parameters
***************************************************************************/
void GModelSpectralTablePars::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone table model parameters
***************************************************************************/
GModelSpectralTablePars* GModelSpectralTablePars::clone(void) const
{
    // Clone table model parameters
    return new GModelSpectralTablePars(*this);
}


/***********************************************************************//**
 * @brief Return pointer to table model parameter
 *
 * @param[in] index Parameter index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Table model parameter index is out of range.
 *
 * Returns a pointer to the table model parameter with the specified
 * @p index.
 ***************************************************************************/
GModelSpectralTablePar* GModelSpectralTablePars::at(const int& index)
{
    // Compile option: raise an exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT, "Parameter index", index, size());
    }
    #endif

    // Return pointer
    return m_pars[index];
}


/***********************************************************************//**
 * @brief Return pointer to table model parameter (const version)
 *
 * @param[in] index Parameter index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Table model parameter index is out of range.
 *
 * Returns a pointer to the table model parameter with the specified
 * @p index.
 ***************************************************************************/
const GModelSpectralTablePar* GModelSpectralTablePars::at(const int& index) const
{
    // Compile option: raise an exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_AT, "Parameter index", index, size());
    }
    #endif

    // Return pointer
    return m_pars[index];
}


/***********************************************************************//**
 * @brief Set table model parameter in container
 *
 * @param[in] index Table model parameter index [0,...,size()-1].
 * @param[in] par Table model parameter.
 * @return Pointer to deep copy of table model parameter.
 *
 * @exception GException::out_of_range
 *            Table model parameter index is out of range.
 * @exception GException::invalid_value
 *            Name of table model parameter exists already in container.
 *
 * Set table model parameter at position @p index in the container. The
 * method will overwrite the table model parameter that existed in the
 * specified slot. The method will store a deep copy of the table model
 * parameter in the container.
 ***************************************************************************/
GModelSpectralTablePar* GModelSpectralTablePars::set(const int&                    index,
                                                     const GModelSpectralTablePar& par)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_SET1, index, size());
    }
    #endif

    // Check if a parameter with specified name does not yet exist
    int inx = get_index(par.par().name());
    if (inx != -1 && inx != index) {
        std::string msg =
            "Attempt to set table model parameter with name \""+
            par.par().name()+"\" in table model parameter container at index "+
            gammalib::str(index)+", but a table model parameter with the "
            "same name exists already at index "+gammalib::str(inx)+
            " in the container. Every table model parameter in the container "
            "needs a unique name.";
        throw GException::invalid_value(G_SET1, msg);
    }

    // Free existing table model parameter only if it differs from current
    // table model parameter. This prevents unintential deallocation of the
    // argument
    if ((m_pars[index] != NULL) && (m_pars[index] != &par)) {
        delete m_pars[index];
    }

    // Assign new table model parameter by cloning
    m_pars[index] = par.clone();

    // Return pointer to table model parameter
    return m_pars[index];
}


/***********************************************************************//**
 * @brief Set table model parameter in container
 *
 * @param[in] name Table model parameter name.
 * @param[in] par Table model parameter.
 * @return Pointer to deep copy of table model parameter.
 *
 * @exception GException::invalid_argument
 *            Table model parameter with @p name not found.
 *
 * Set table model parameter at the position of the parameter with the
 * specified @p name in the container. The method will overwrite the table
 * model parameter with the specified @p name. The method will store a deep
 * copy of the table model parameter in the container.
 ***************************************************************************/
GModelSpectralTablePar* GModelSpectralTablePars::set(const std::string&            name,
                                                     const GModelSpectralTablePar& par)
{
    // Get parameter index
    int index = get_index(name);

    // Throw exception if model name was not found
    if (index == -1) {
        std::string msg = "Table model parameter \""+name+"\" not found. "
                          "Please specify a valid parameter name.";
        throw GException::invalid_argument(G_SET2, msg);
    }

    // Set parameter and return pointer to table model parameter
    return (set(index, par));
}


/***********************************************************************//**
 * @brief Append table model parameter to container
 *
 * @param[in] par Table model parameter.
 * @return Pointer to deep copy of table model parameter.
 *
 * @exception GException::invalid_argument
 *            Name of table model parameter exists already in container.
 *
 * Appends table model parameter to the container by making a deep copy of
 * the table model parameter and storing its pointer.
 ***************************************************************************/
GModelSpectralTablePar* GModelSpectralTablePars::append(const GModelSpectralTablePar& par)
{
    // Check if a model with specified name does not yet exist
    int inx = get_index(par.par().name());
    if (inx != -1) {
        std::string msg = 
            "Attempt to append table model parameter with name \""+
            par.par().name()+"\" to table model parameter container, but a "
            "table model parameter with the same name exists already at "
            "index "+gammalib::str(inx)+" in the container. Every table "
            "model parameter in the container needs a unique name.";
        throw GException::invalid_argument(G_APPEND, msg);
    }

    // Create deep copy of table model parameter
    GModelSpectralTablePar* ptr = par.clone();

    // Append deep copy of table model parameter
    m_pars.push_back(ptr);

    // Return pointer to table model parameter
    return ptr;
}


/***********************************************************************//**
 * @brief Insert table model parameter into container
 *
 * @param[in] index Table model parameter index [0,...,size()-1].
 * @param[in] par Table model parameter.
 * @return Pointer to deep copy of table model parameter.
 *
 * @exception GException::out_of_range
 *            Model index is out of range.
 * @exception GException::invalid_value
 *            Name of model exists already in container.
 *
 * Inserts a @p model into the container before the model with the specified
 * @p index.
 ***************************************************************************/
GModelSpectralTablePar* GModelSpectralTablePars::insert(const int&                    index,
                                                        const GModelSpectralTablePar& par)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (is_empty()) {
        if (index > 0) {
            throw GException::out_of_range(G_INSERT1, index, size());
        }
    }
    else {
        if (index < 0 || index >= size()) {
            throw GException::out_of_range(G_INSERT1, index, size());
        }
    }
    #endif

    // Check if a model with specified name does not yet exist
    int inx = get_index(par.par().name());
    if (inx != -1) {
        std::string msg =
            "Attempt to insert table model parameter with name \""+
            par.par().name()+"\" in table model parameter container before "
            "index "+gammalib::str(index)+", but table model parameter with "
            "the same name exists already at index "+
            gammalib::str(inx)+" in the container. Every model in the table "
            "model parameter container needs a unique name.";
        throw GException::invalid_value(G_INSERT1, msg);
    }

    // Create deep copy of table model parameter
    GModelSpectralTablePar* ptr = par.clone();

    // Inserts deep copy of table model parameter
    m_pars.insert(m_pars.begin()+index, ptr);

    // Return pointer to table model parameter
    return ptr;
}


/***********************************************************************//**
 * @brief Insert table model parameter into container
 *
 * @param[in] name Table model parameter name.
 * @param[in] par Table model parameter.
 * @return Pointer to deep copy of table model parameter.
 *
 * @exception GException::invalid_argument
 *            Table model parameter with @p name not found.
 *
 * Inserts a table model parameter into the container before the table model
 * parameter with the specified @p name.
 ***************************************************************************/
GModelSpectralTablePar* GModelSpectralTablePars::insert(const std::string&            name,
                                                        const GModelSpectralTablePar& par)
{
    // Get parameter index
    int index = get_index(name);

    // Throw exception if model name was not found
    if (index == -1) {
        std::string msg = "Table model parameter \""+name+"\" not found. "
                          "Please specify a valid parameter name.";
        throw GException::invalid_argument(G_INSERT2, msg);
    }

    // Insert table model parameter and return pointer to it
    return (insert(index, par));
}


/***********************************************************************//**
 * @brief Remove table model parameter from container
 *
 * @param[in] index Table model parameter index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Table model parameter index is out of range.
 *
 * Remove table model parameter of specified @p index from container.
 ***************************************************************************/
void GModelSpectralTablePars::remove(const int& index)
{
    // Compile option: raise exception if index is out of range
    #if defined(G_RANGE_CHECK)
    if (index < 0 || index >= size()) {
        throw GException::out_of_range(G_REMOVE1, index, size());
    }
    #endif

    // Delete model
    delete m_pars[index];

    // Erase model component from container
    m_pars.erase(m_pars.begin() + index);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Remove table model parameter from container
 *
 * @param[in] name Table model parameter name.
 *
 * @exception GException::invalid_argument
 *            Table model parameter with @p name not found.
 *
 * Remove table model parameter with @p name from container.
 ***************************************************************************/
void GModelSpectralTablePars::remove(const std::string& name)
{
    // Get parameter index
    int index = get_index(name);

    // Throw exception if model name was not found
    if (index == -1) {
        std::string msg = "Table model parameter \""+name+"\" not found. "
                          "Please specify a valid parameter name.";
        throw GException::invalid_argument(G_REMOVE2, msg);
    }

    // Remove model
    remove(index);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Append table model parameter container
 *
 * @param[in] pars Table model parameter container.
 *
 * Append table model parameter container to the container.
 ***************************************************************************/
void GModelSpectralTablePars::extend(const GModelSpectralTablePars& pars)
{
    // Do nothing if table model parameter container is empty
    if (!pars.is_empty()) {

        // Get size. Note that we extract the size first to avoid an
        // endless loop that arises when a container is appended to
        // itself.
        int num = pars.size();

        // Reserve enough space
        reserve(size() + num);

        // Loop over all table model parameters and append pointers to deep
        // copies
        for (int i = 0; i < num; ++i) {

            // Check if table model parameter name does not yet exist
            int inx = get_index(pars[i]->par().name());
            if (inx != -1) {
                std::string msg =
                    "Attempt to append table model parameter with name \""+
                    pars[i]->par().name()+"\" to table model parameter "
                    "container, but a table model parameter with the same "
                    "name exists already at index "+gammalib::str(inx)+" "
                    "in the container. Every table model parameter in the "
                    "container needs a unique name.";
                throw GException::invalid_value(G_EXTEND, msg);
            }

            // Append model to container
            m_pars.push_back(pars[i]->clone());

        } // endfor: looped over all models

    } // endif: model container was not empty
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Signals if table model parameter exists
 *
 * @param[in] name Table model parameter name.
 * @return True if table model parameter with specified @p name exists.
 *
 * Searches all table model parameters for a match with the specified
 * @p name. If the specified name has been found, true is returned.
 ***************************************************************************/
bool GModelSpectralTablePars::contains(const std::string& name) const
{
    // Get model index
    int index = get_index(name);

    // Return
    return (index != -1);
}


/***********************************************************************//**
 * @brief Print table model parameters
 *
 * @param[in] chatter Chattiness.
 * @return String containing file table model parameters information.
 ***************************************************************************/
std::string GModelSpectralTablePars::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GModelSpectralTablePars ===");

        // Append information
        result.append("\n"+gammalib::parformat("Number of parameters"));
        result.append(gammalib::str(size()));
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
void GModelSpectralTablePars::init_members(void)
{
    // Initialise members
    m_pars.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] pars Table model parameters.
 ***************************************************************************/
void GModelSpectralTablePars::copy_members(const GModelSpectralTablePars& pars)
{
    // Copy parameters
    m_pars.clear();
    for (int i = 0; i < pars.m_pars.size(); ++i) {
        m_pars.push_back((pars.m_pars[i]->clone()));
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GModelSpectralTablePars::free_members(void)
{
    // Free models
    for (int i = 0; i < m_pars.size(); ++i) {
        delete m_pars[i];
        m_pars[i] = NULL;
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
int GModelSpectralTablePars::get_index(const std::string& name) const
{
    // Initialise index
    int index = -1;

    // Search model with specified name
    for (int i = 0; i < size(); ++i) {
        if (m_pars[i]->par().name() == name) {
            index = i;
            break;
        }
    }

    // Return index
    return index;
}
