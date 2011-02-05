/***************************************************************************
 *                   GPars.hpp - Application parameters                    *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GPars.hpp
 * @brief Application parameters class definition
 * @author Jurgen Knodlseder
 */

#ifndef GPARS_HPP
#define GPARS_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include <iostream>
#include "GPar.hpp"
#include "GLog.hpp"


/***********************************************************************//**
 * @class GPars
 *
 * @brief Application parameters interface defintion.
 ***************************************************************************/
class GPars {

    // Friend classes
    friend class GApplication;

    // I/O friends
    friend std::ostream& operator<<(std::ostream& os, const GPars& pars);
    friend GLog&         operator<<(GLog&        log, const GPars& pars);

public:
    // Constructors and destructors
    GPars(void);
    GPars(const std::string& filename);
    GPars(const std::string& filename, const std::vector<std::string>& args);
    GPars(const GPars& pars);
    ~GPars(void);
 
    // Operators
    GPars& operator= (const GPars& pars);

    // Methods
    void        clear(void);
    int         size(void) const { return m_pars.size(); }
    void        load(const std::string& filename);
    void        load(const std::string& filename, const std::vector<std::string>& args);
    void        save(const std::string& filename);
    GPar*       par(const std::string& name);
    const GPar* par(const std::string& name) const;
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
