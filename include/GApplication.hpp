/***************************************************************************
 *             GApplication.hpp - GammaLib application base class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GApplication.hpp
 * @brief GammaLib application base class
 * @author Jurgen Knodlseder
 */

#ifndef GAPPLICATION_HPP
#define GAPPLICATION_HPP

/* __ Includes ___________________________________________________________ */
#include <time.h>
#include <vector>
#include <string>
#include <iostream>
#include "GLog.hpp"
#include "GPars.hpp"


/***********************************************************************//**
 * @class GApplication
 *
 * @brief GammaLib application interface defintion.
 ***************************************************************************/
class GApplication {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GApplication& app);

public:
    // Constructors and destructors
    GApplication(void);
    GApplication(const std::string& name, const std::string& version);
    GApplication(const std::string& name, const std::string& version,
                 int argc, char* argv[]);
    GApplication(const GApplication& app);
    ~GApplication(void);
 
    // Operators
    GApplication& operator= (const GApplication& app);

    // Methods
    std::string name(void) const;
    std::string version(void) const;
    double      telapse(void) const;
  
protected:
    // Protected methods
    void        init_members(void);
    void        copy_members(const GApplication& app);
    void        free_members(void);
    std::string parfilename(void) const;

    // Protected data members
    std::string              m_name;       //!< Application name
    std::string              m_version;    //!< Application version
    std::vector<std::string> m_args;       //!< Command line arguments
    time_t                   m_tstart;     //!< Start time of execution
    GLog                     m_logger;     //!< Application logger
    GPars                    m_pars;       //!< Application parameters
};

#endif /* GAPPLICATION_HPP */
