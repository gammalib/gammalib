/***************************************************************************
 *                   GPars.hpp - Application parameters                    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2012 by Juergen Knoedlseder                         *
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
#include <iostream>
#include "GBase.hpp"
#include "GPar.hpp"
#include "GLog.hpp"


/***********************************************************************//**
 * @class GPars
 *
 * @brief Application parameter container class
 *
 * This class holds a collection of application parameters.
 ***************************************************************************/
class GPars : public GBase {

    // Friend classes
    friend class GApplication;

    // I/O friends
    friend std::ostream& operator<<(std::ostream& os, const GPars& pars);
    friend GLog&         operator<<(GLog&        log, const GPars& pars);

public:
    // Constructors and destructors
    GPars(void);
    explicit GPars(const std::string& filename);
    explicit GPars(const std::string& filename, const std::vector<std::string>& args);
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
    int         size(void) const { return m_pars.size(); }
    void        append(const GPar& par);
    void        append_standard(void);
    void        load(const std::string& filename);
    void        load(const std::string& filename, const std::vector<std::string>& args);
    void        save(const std::string& filename);
    bool        haspar(const std::string& name) const;
    std::string print(void) const;
  
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

    // Protected data members
    std::vector<std::string> m_parfile;   //!< Parameter file lines
    std::vector<GPar>        m_pars;      //!< Parameters
    std::vector<int>         m_line;      //!< Line number of parameter
    std::vector<size_t>      m_vstart;    //!< Column of value start
    std::vector<size_t>      m_vstop;     //!< Column of value stop
    std::string              m_mode;      //!< Effective mode
};

#endif /* GPARS_HPP */
