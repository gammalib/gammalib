/***************************************************************************
 *                   GPar.hpp - Application parameter class                *
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
 * @file GPar.hpp
 * @brief Application parameter class definition
 * @author Juergen Knoedlseder
 */

#ifndef GPAR_HPP
#define GPAR_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"


/***********************************************************************//**
 * @class GPar
 *
 * @brief Application parameter class interface defintion.
 ***************************************************************************/
class GPar : public GBase {

    // Friend classes
    friend class GPars;
    friend class GApplication;

public:
    // Constructors and destructors
    GPar(void);
    explicit GPar(const std::string& name, const std::string& type,
                  const std::string& mode, const std::string& value,
                  const std::string& min, const std::string& max, 
                  const std::string& prompt);
    GPar(const GPar& par);
    virtual ~GPar(void);
 
    // Operators
    GPar& operator= (const GPar& par);

    // Methods
    void        clear(void);
    GPar*       clone(void) const;
    void        type(const std::string& type);
    void        mode(const std::string& mode);
    void        value(const std::string& value);
    void        string(const std::string& value);
    void        filename(const std::string& value);
    void        boolean(const bool& value);
    void        integer(const int& value);
    void        real(const double& value);
    std::string name(void) const { return m_name; }
    std::string type(void) const { return m_type; }
    std::string mode(void) const { return m_mode; }
    std::string value(void);
    std::string string(void);
    std::string filename(void);
    bool        boolean(void);
    int         integer(void);
    double      real(void);
    std::string min(void) const { return m_min; }
    std::string max(void) const { return m_max; }
    std::string prompt(void) const { return m_prompt; }
    bool        islearn(void) const;
    bool        isquery(void) const;
    bool        isfilename(void) const;
    std::string print(void) const;
  
protected:
    // Protected methods
    void        init_members(void);
    void        copy_members(const GPar& par);
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

    // Protected data members
    bool        m_update;  //!< Signal value updating
    std::string m_name;    //!< Parameter name
    std::string m_type;    //!< Parameter type
    std::string m_mode;    //!< Parameter mode
    std::string m_value;   //!< Parameter value
    std::string m_min;     //!< Parameter minimum
    std::string m_max;     //!< Parameter maximum
    std::string m_prompt;  //!< Parameter prompt
};

#endif /* GPAR_HPP */
