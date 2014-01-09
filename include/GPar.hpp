/***************************************************************************
 *                   GApplicationPar.hpp - Application parameter class                *
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
 * @file GApplicationPar.hpp
 * @brief Application parameter class definition
 * @author Juergen Knoedlseder
 */

#ifndef GPAR_HPP
#define GPAR_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"


/***********************************************************************//**
 * @class GApplicationPar
 *
 * @brief Application parameter class
 ***************************************************************************/
class GApplicationPar : public GBase {

    // Friend classes
    friend class GApplicationPars;
    friend class GApplication;

public:
    // Constructors and destructors
    GApplicationPar(void);
    explicit GApplicationPar(const std::string& name, const std::string& type,
                  const std::string& mode, const std::string& value,
                  const std::string& min, const std::string& max, 
                  const std::string& prompt);
    GApplicationPar(const GApplicationPar& par);
    virtual ~GApplicationPar(void);
 
    // Operators
    GApplicationPar& operator=(const GApplicationPar& par);

    // Methods
    void               clear(void);
    GApplicationPar*              clone(void) const;
    void               type(const std::string& type);
    void               mode(const std::string& mode);
    void               value(const std::string& value);
    void               string(const std::string& value);
    void               filename(const std::string& value);
    void               boolean(const bool& value);
    void               integer(const int& value);
    void               real(const double& value);
    const std::string& name(void) const;
    const std::string& type(void) const;
    const std::string& mode(void) const;
    std::string        value(void);
    std::string        string(void);
    std::string        filename(void);
    bool               boolean(void);
    int                integer(void);
    double             real(void);
    const std::string& min(void) const;
    const std::string& max(void) const;
    const std::string& prompt(void) const;
    bool               is_learn(void) const;
    bool               is_query(void) const;
    bool               is_filename(void) const;
    bool               is_valid(void);
    bool               is_undefined(void);
    bool               is_notanumber(void);
    std::string        print(const GChatter& chatter = NORMAL) const;
  
protected:
    // Protected enumerators
    enum Status {
        ST_VALID,
        ST_UNDEFINED,
        ST_NAN,
        ST_UNDERFLOW,
        ST_OVERFLOW
    };

    // Protected methods
    void        init_members(void);
    void        copy_members(const GApplicationPar& par);
    void        free_members(void);
    void        check_type(const std::string& type) const;
    void        check_mode(const std::string& mode) const;
    void        check_value(const std::string& value) const;
    void        check_value_bool(const std::string& value) const;
    void        check_value_int(const std::string& value) const;
    void        check_value_real(const std::string& value) const;
    void        check_value_string(const std::string& value) const;
    void        check_value_filename(const std::string& value) const;
    void        set_value(const std::string& value);
    void        query(void);
    void        stop_query(void);
    std::string par_type_string(const std::string& type) const;
    std::string par_status_string(void) const;

    // Protected data members
    bool        m_update;  //!< Signal value updating
    std::string m_name;    //!< Parameter name
    std::string m_type;    //!< Parameter type
    std::string m_mode;    //!< Parameter mode
    std::string m_value;   //!< Parameter value
    std::string m_min;     //!< Parameter minimum
    std::string m_max;     //!< Parameter maximum
    std::string m_prompt;  //!< Parameter prompt
    Status      m_status;  //!< Parameter status
};


/***********************************************************************//**
 * @brief Returns parameter name
 *
 * @return Parameter name
 ***************************************************************************/
inline
const std::string& GApplicationPar::name(void) const
{
    return m_name;
}


/***********************************************************************//**
 * @brief Returns parameter type
 *
 * @return Parameter type
 ***************************************************************************/
inline
const std::string& GApplicationPar::type(void) const
{
    return m_type;
}


/***********************************************************************//**
 * @brief Returns parameter mode
 *
 * @return Parameter mode
 ***************************************************************************/
inline
const std::string& GApplicationPar::mode(void) const
{
    return m_mode;
}


/***********************************************************************//**
 * @brief Returns parameter minimum
 *
 * @return Parameter minimum
 ***************************************************************************/
inline
const std::string& GApplicationPar::min(void) const
{
    return m_min;
}


/***********************************************************************//**
 * @brief Returns parameter maximum
 *
 * @return Parameter maximum
 ***************************************************************************/
inline
const std::string& GApplicationPar::max(void) const
{
    return m_max;
}


/***********************************************************************//**
 * @brief Returns parameter prompt
 *
 * @return Parameter prompt
 ***************************************************************************/
inline
const std::string& GApplicationPar::prompt(void) const
{
    return m_prompt;
}

#endif /* GPAR_HPP */
