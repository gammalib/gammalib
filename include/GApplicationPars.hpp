/***************************************************************************
 *               GApplicationPars.hpp - Application parameters             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2018 by Juergen Knoedlseder                         *
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
 * @file GApplicationPars.hpp
 * @brief Application parameter container class definition
 * @author Juergen Knoedlseder
 */

#ifndef GAPPLICATIONPARS_HPP
#define GAPPLICATIONPARS_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include "GContainer.hpp"
#include "GApplicationPar.hpp"

/* __ Forward declarations _______________________________________________ */
class GFilename;


/***********************************************************************//**
 * @class GApplicationPars
 *
 * @brief Application parameter container class
 *
 * This class holds a collection of application parameters.
 ***************************************************************************/
class GApplicationPars : public GContainer {

    // Friend classes
    friend class GApplication;

public:
    // Constructors and destructors
    GApplicationPars(void);
    explicit GApplicationPars(const GFilename& filename);
    GApplicationPars(const GFilename&                filename,
                     const std::vector<std::string>& args);
    GApplicationPars(const GApplicationPars& pars);
    virtual ~GApplicationPars(void);
 
    // Operators
    GApplicationPars&      operator=(const GApplicationPars& pars);
    GApplicationPar&       operator[](const int& index);
    const GApplicationPar& operator[](const int& index) const;
    GApplicationPar&       operator[](const std::string& name);
    const GApplicationPar& operator[](const std::string& name) const;

    // Methods
    void                     clear(void);
    GApplicationPars*        clone(void) const;
    std::string              classname(void) const;
    GApplicationPar&         at(const int& index);
    const GApplicationPar&   at(const int& index) const;
    int                      size(void) const;
    bool                     is_empty(void) const;
    GApplicationPar&         append(const GApplicationPar& par);
    void                     append_standard(void);
    GApplicationPar&         insert(const int& index, const GApplicationPar& par);
    GApplicationPar&         insert(const std::string&     name,
                                    const GApplicationPar& par);
    void                     remove(const int& index);
    void                     remove(const std::string& name);
    void                     reserve(const int& num);
    void                     extend(const GApplicationPars& pars);
    bool                     contains(const std::string& name) const;
    void                     syspfiles(const std::string& syspfiles);
    const std::string&       syspfiles(void) const;
    void                     load(const GFilename& filename);
    void                     load(const GFilename&                filename,
                                  const std::vector<std::string>& args);
    void                     save(const GFilename& filename);
    void                     pickle(const std::vector<std::string>& string);
    std::vector<std::string> pickle(void) const;
    std::string              print(const GChatter& chatter = NORMAL) const;
  
protected:
    // Protected methods
    void        init_members(void);
    void        copy_members(const GApplicationPars& pars);
    void        free_members(void);
    std::string inpath(const std::string& filename) const;
    std::string pfiles_path(const std::string& filename) const;
    std::string syspfiles_path(const std::string& filename) const;
    std::string outpath(const std::string& filename) const;
    void        read(const std::string& filename);
    void        write(const std::string& filename) const;
    void        parse(void);
    void        update(void);
    void        synchronise(const std::string& filename);
    int         get_index(const std::string& name) const;
    std::string parline(GApplicationPar& par, size_t* start, size_t* stop) const;

    // Protected data members
    std::vector<std::string>     m_parfile;   //!< Parameter file lines
    std::vector<GApplicationPar> m_pars;      //!< Parameters
    std::vector<int>             m_line;      //!< Line number of parameter
    std::vector<size_t>          m_vstart;    //!< Column of value start
    std::vector<size_t>          m_vstop;     //!< Column of value stop
    std::string                  m_mode;      //!< Effective mode
    std::string                  m_syspfiles; //!< Optional location of syspfiles
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GApplicationPars").
 ***************************************************************************/
inline
std::string GApplicationPars::classname(void) const
{
    return ("GApplicationPars");
}


/***********************************************************************//**
 * @brief Returns reference to parameter
 *
 * @param[in] index Parameter index [0,...,size()-1].
 *
 * Returns a reference to the parameter with the specified @p index.
 ***************************************************************************/
inline
GApplicationPar& GApplicationPars::operator[](const int& index)
{
    return (m_pars[index]);
}


/***********************************************************************//**
 * @brief Returns reference to parameter (const version)
 *
 * @param[in] index Parameter index [0,...,size()-1].
 *
 * Returns a reference to the parameter with the specified @p index.
 ***************************************************************************/
inline
const GApplicationPar& GApplicationPars::operator[](const int& index) const
{
    return (m_pars[index]);
}


/***********************************************************************//**
 * @brief Return number of parameters in container
 *
 * @return Number of parameters in container.
 *
 * Returns the number of parameters in the parameter container.
 ***************************************************************************/
inline
int GApplicationPars::size(void) const
{
    return (int)m_pars.size();
}


/***********************************************************************//**
 * @brief Signals if there are no parameters in container
 *
 * @return True if container is empty, false otherwise.
 *
 * Signals if the parameter container does not contain any parameter.
 ***************************************************************************/
inline
bool GApplicationPars::is_empty(void) const
{
    return (m_pars.empty());
}


/***********************************************************************//**
 * @brief Reserves space for parameters in container
 *
 * @param[in] num Number of parameters.
 *
 * Reserves space for @p num parameters in the container.
 ***************************************************************************/
inline
void GApplicationPars::reserve(const int& num)
{
    m_pars.reserve(num);
    return;
}


/***********************************************************************//**
 * @brief Set path to system pfiles
 *
 * @param[in] syspfiles Path to system pfiles.
 *
 * Sets the path to system parameter files.
 ***************************************************************************/
inline
void GApplicationPars::syspfiles(const std::string& syspfiles)
{
    m_syspfiles = syspfiles;
    return;
}


/***********************************************************************//**
 * @brief Return path to system pfiles
 *
 * @return Path to system pfiles.
 *
 * Returns path to system parameter files.
 ***************************************************************************/
inline
const std::string& GApplicationPars::syspfiles(void) const
{
    return (m_syspfiles);
}

#endif /* GAPPLICATIONPARS_HPP */
