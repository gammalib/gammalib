/***************************************************************************
 *                   GPars.hpp - Application parameters                    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2013 by Juergen Knoedlseder                         *
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
 * @file GPars.hpp
 * @brief Application parameter container class definition
 * @author Juergen Knoedlseder
 */

#ifndef GPARS_HPP
#define GPARS_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include "GContainer.hpp"
#include "GPar.hpp"


/***********************************************************************//**
 * @class GPars
 *
 * @brief Application parameter container class
 *
 * This class holds a collection of application parameters.
 ***************************************************************************/
class GPars : public GContainer {

    // Friend classes
    friend class GApplication;

public:
    // Constructors and destructors
    GPars(void);
    explicit GPars(const std::string& filename);
    explicit GPars(const std::string& filename,
                   const std::vector<std::string>& args);
    GPars(const GPars& pars);
    virtual ~GPars(void);
 
    // Operators
    GPars&      operator=(const GPars& pars);
    GPar&       operator[](const int& index);
    const GPar& operator[](const int& index) const;
    GPar&       operator[](const std::string& name);
    const GPar& operator[](const std::string& name) const;

    // Methods
    void        clear(void);
    GPars*      clone(void) const;
    GPar&       at(const int& index);
    const GPar& at(const int& index) const;
    int         size(void) const;
    bool        isempty(void) const;
    GPar&       append(const GPar& par);
    void        append_standard(void);
    GPar&       insert(const int& index, const GPar& par);
    GPar&       insert(const std::string& name, const GPar& par);
    void        remove(const int& index);
    void        remove(const std::string& name);
    void        reserve(const int& num);
    void        extend(const GPars& pars);
    bool        contains(const std::string& name) const;
    void        load(const std::string& filename);
    void        load(const std::string& filename,
                     const std::vector<std::string>& args);
    void        save(const std::string& filename);
    std::string print(const GChatter& chatter = NORMAL) const;
  
protected:
    // Protected methods
    void        init_members(void);
    void        copy_members(const GPars& pars);
    void        free_members(void);
    std::string inpath(const std::string& filename) const;
    std::string outpath(const std::string& filename) const;
    void        read(const std::string& filename);
    void        write(const std::string& filename) const;
    void        parse(void);
    void        update(void);
    int         get_index(const std::string& name) const;
    std::string parline(GPar& par, size_t* start, size_t* stop) const;

    // Protected data members
    std::vector<std::string> m_parfile;   //!< Parameter file lines
    std::vector<GPar>        m_pars;      //!< Parameters
    std::vector<int>         m_line;      //!< Line number of parameter
    std::vector<size_t>      m_vstart;    //!< Column of value start
    std::vector<size_t>      m_vstop;     //!< Column of value stop
    std::string              m_mode;      //!< Effective mode
};


/***********************************************************************//**
 * @brief Returns reference to parameter
 *
 * @param[in] index Parameter index [0,...,size()-1].
 *
 * Returns a reference to the parameter with the specified @p index.
 ***************************************************************************/
inline
GPar& GPars::operator[](const int& index)
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
const GPar& GPars::operator[](const int& index) const
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
int GPars::size(void) const
{
    return (m_pars.size());
}


/***********************************************************************//**
 * @brief Signals if there are no parameters in container
 *
 * @return True if container is empty, false otherwise.
 *
 * Signals if the parameter container does not contain any parameter.
 ***************************************************************************/
inline
bool GPars::isempty(void) const
{
    return (m_pars.empty());
}


/***********************************************************************//**
 * @brief Reserves space for parameters in container
 *
 * @param[in] num Number of parameter
 *
 * Reserves space for @p num parameters in the container.
 ***************************************************************************/
inline
void GPars::reserve(const int& num)
{
    m_pars.reserve(num);
    return;
}

#endif /* GPARS_HPP */
