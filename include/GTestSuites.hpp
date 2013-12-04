/***************************************************************************
 *              GTestSuites.hpp - Test suite container class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2013 by Jean-Baptiste Cayrou                        *
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
 * @file GTestSuites.hpp
 * @brief Test suite container class interface defintion
 * @author Jean-Baptiste Cayrou
 */

#ifndef GTESTSUITES_HPP
#define GTESTSUITES_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GContainer.hpp"
#include "GLog.hpp"
#include "GXml.hpp"
#include "GTestSuite.hpp"


/***********************************************************************//**
 * @class GTestSuites
 *
 * @brief Test suite container class
 *
 * This is a container class for test suites. All test suites will be
 * executed by invoking the run() method, and test results are saved in
 * JUnit format to an XML file using the save() method.
 *
 * @todo Provide detailed description
 ***************************************************************************/
class GTestSuites : public GContainer {
    
public:
    // Constructors and destructors
    GTestSuites(void);
    GTestSuites(const GTestSuites& suites);
    GTestSuites(const std::string& name);
    virtual ~GTestSuites(void);
    
    // Operators
    GTestSuites&      operator=(const GTestSuites& suites);
    GTestSuite*       operator[](const int& index);
    const GTestSuite* operator[](const int& index) const;
        
    // Methods
    void               clear(void);
    GTestSuites*       clone(void) const;
    GTestSuite*        at(const int& index);
    const GTestSuite*  at(const int& index) const;
    int                size(void) const;
    bool               isempty(void) const;
    GTestSuite*        set(const int& index, const GTestSuite& suite);
    GTestSuite*        append(const GTestSuite& suite);
    GTestSuite*        insert(const int& index, const GTestSuite& suite);
    void               remove(const int& index);
    void               reserve(const int& num);
    void               extend(const GTestSuites& suites);
    const std::string& name(void) const;
    void               name(const std::string& name);
    void               cout(const bool& flag);
    int                errors(void) const;
    int                failures(void) const;
    int                tests(void) const;
    const time_t&      timestamp(void) const;
    bool               run(void);
    void               save(const std::string& filename) const;
    std::string        print(const GChatter& chatter = NORMAL) const;
    
protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GTestSuites& suites);
    void free_members(void);
    void write(GXml& xml) const;
        
    // Protected members    
    std::string              m_name;       //!< Name of container
    std::vector<GTestSuite*> m_testsuites; //!< Vector of test suites
    time_t                   m_timestamp;  //!< Timestamp
    GLog                     m_log;        //!< Log
};


/***********************************************************************//**
 * @brief Returns pointer to test suite
 *
 * @param[in] index Test suite index [0,...,size()-1].
 * @return Pointer to test suite.
 *
 * Returns pointer to test suite with @p index. No range checking is done
 * for the index.
 ***************************************************************************/
inline
GTestSuite* GTestSuites::operator[](const int& index)
{
    return (m_testsuites[index]);
}


/***********************************************************************//**
 * @brief Returns pointer to test suite (const version)
 *
 * @param[in] index Test suite index [0,...,size()-1].
 * @return Pointer to test suite.
 *
 * Returns pointer to test suite with @p index. No range checking is done
 * for the index.
 ***************************************************************************/
inline
const GTestSuite* GTestSuites::operator[](const int& index) const
{
    return (m_testsuites[index]);
}


/***********************************************************************//**
 * @brief Return number of test suites in container
 *
 * @return Number of test suites in container.
 *
 * Returns the number of test suites in the container.
 ***************************************************************************/
inline
int GTestSuites::size(void) const
{
    return m_testsuites.size();
}


/***********************************************************************//**
 * @brief Signals if there are no test suites in the container
 *
 * @return True if test suite container is empty, false otherwise.
 *
 * Signals if the container does not contain any test suite.
 ***************************************************************************/
inline
bool GTestSuites::isempty(void) const
{
    return (m_testsuites.empty());
}


/***********************************************************************//**
 * @brief Reserves space for test suites in container
 *
 * @param[in] num Number of test suites.
 *
 * Reserves space for @p num test suites in the container.
 ***************************************************************************/
inline
void GTestSuites::reserve(const int& num)
{
    m_testsuites.reserve(num);
    return;
}


/***********************************************************************//**
 * @brief Return test suite container name
 *
 * @return Test suite container name
 ***************************************************************************/
inline
const std::string& GTestSuites::name(void) const
{
    return m_name;
}


/***********************************************************************//**
 * @brief Set test suite container name
 *
 * @param[in] name Test suite container name.
 ***************************************************************************/
inline
void GTestSuites::name(const std::string& name)
{
    m_name = name;
    return;
}


/***********************************************************************//**
 * @brief Enables/disables logging into standard output stream
 *
 * @param[in] flag Enable/disable logging (true/false).
 *
 * Enables or disables logging into the standard output stream.
 ***************************************************************************/
inline
void GTestSuites::cout(const bool& flag)
{
    m_log.cout(flag);
    return;
}


/***********************************************************************//**
 * @brief Return the timestamp
 *
 * The time step is set at the moment of construction of the test suites
 * container.
 ***************************************************************************/
inline
const time_t& GTestSuites::timestamp(void) const
{
    return m_timestamp;
}

#endif /* GTESTSUITES_HPP */
