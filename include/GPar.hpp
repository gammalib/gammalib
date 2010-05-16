/***************************************************************************
 *                GPar.hpp - Application parameter base class              *
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
 * @file GPar.hpp
 * @brief Application parameter base class definition
 * @author Jurgen Knodlseder
 */

#ifndef GPAR_HPP
#define GPAR_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <iostream>


/***********************************************************************//**
 * @class GPar
 *
 * @brief Application parameter base class interface defintion.
 ***************************************************************************/
class GPar {

    // Friend classes
    friend class GPars;

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GPar& par);

public:
    // Constructors and destructors
    GPar(void);
    GPar(const std::string& name, const std::string& type,
         const std::string& mode, const std::string& value,
         const std::string& min, const std::string& max, 
         const std::string& desc);
    GPar(const GPar& par);
    ~GPar(void);
 
    // Operators
    GPar& operator= (const GPar& par);

    // Methods
  
protected:
    // Protected methods
    void  init_members(void);
    void  copy_members(const GPar& par);
    void  free_members(void);
    GPar* clone(void) const;

    // Protected data members
    std::string m_name;    //!< Parameter name
    std::string m_type;    //!< Parameter type
    std::string m_mode;    //!< Parameter mode
    std::string m_value;   //!< Parameter value
    std::string m_min;     //!< Parameter minimum
    std::string m_max;     //!< Parameter maximum
    std::string m_desc;    //!< Parameter descriptor
};

#endif /* GPAR_HPP */
