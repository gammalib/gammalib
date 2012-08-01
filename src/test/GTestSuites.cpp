/***************************************************************************
 *         GTestSuites.cpp  - Test Suites class for GammaLib               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012 by Jean-Baptiste Cayrou                             *
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
 * @file GTestSuites.cpp
 * @brief Test Suites class implementation
 * @author Jean-Baptiste Cayrou
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include<ctime>
#include "GTestSuites.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_OP_ACCESS                         "GTestSuites::operator[](int&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GTestSuites::GTestSuites(void)
{
    // Initialise members
    init_members();

    //Return
    return;
}

/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] testsuites Test Suites container.
 ***************************************************************************/
GTestSuites::GTestSuites(const GTestSuites& testsuites)
{
    // Initialise members
    init_members();
    
    // Copy members
    copy_members(testsuites);
    
    //Return
    return;
}

/***********************************************************************//**
 * @brief Name constructor
 *
 * @param[in] name Name of the Test Suites
 ***************************************************************************/
GTestSuites::GTestSuites(const std::string& name)
{
    //Initialise members
    init_members();
    
    //Set name
    m_name=name;

    //Return
    return;
}

/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GTestSuites::~GTestSuites()
{
    //Free members
    free_members();

    //Return
    return;
}

/*==========================================================================
 =                                                                         =
 =                               Operators                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] testsuites Test Suites container.
 ***************************************************************************/
GTestSuites&  GTestSuites::operator= (const GTestSuites& testsuites)
{
    // Execute only if object is not identical
    if (this != &testsuites) {
        
        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(testsuites);

    } // endif: object was not identical

    // Return
    return *this;
}

/***********************************************************************//**
 * @brief Returns reference to test suite
 *
 * @param[in] index Test suites index [0,...,test_suites()-1].
 *
 * @exception GException::out_of_range
 *            Test suites index is out of range.
 ***************************************************************************/
GTestSuite& GTestSuites::operator[](const int& index)
{
    // If index is outside boundary then throw an error
    if (index < 0 || index >= test_suites()) {
        throw GException::out_of_range(G_OP_ACCESS, index, 0, test_suites()-1);
    }

    // Return reference
    return *(m_testsuites[index]);
}

/***********************************************************************//**
 * @brief Returns reference to test suite
 *
 * @param[in] index Test Suites index [0,...,size()-1].
 *
 * @exception GException::out_of_range
 *            Test Suites index is out of range.
 ***************************************************************************/
const GTestSuite& GTestSuites::operator[](const int& index) const
{
    // If index is outside boundary then throw an error
    if (index < 0 || index >= test_suites()) {
        throw GException::out_of_range(G_OP_ACCESS, index, 0, test_suites()-1);
    }

    // Return reference
    return *(m_testsuites[index]);
}

/*==========================================================================
 =                                                                         =
 =                            Public methods                               =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Return Test Suites name
 ***************************************************************************/
std::string GTestSuites::name() const
{
    //return Test Suites name
    return m_name;
}

/***********************************************************************//**
 * @brief Set Test Suites name
 * @param[in] name Parameter name.
 ***************************************************************************/
void GTestSuites::name(const std::string& name)
{
    //Set name
    m_name=name;
    
    //Return
    return;
}

/***********************************************************************//**
 * @brief Append test suite to container
 *
 * @param[in] testsuite Test suite.
 *
 * This method appends one test suite to the container.
 ***************************************************************************/
void GTestSuites::append(GTestSuite& testsuite)
{
    //Add testsuite to container
    m_testsuites.push_back(&testsuite);

    //Return
    return;
}

/***********************************************************************//**
 * @brief Return the number of test suite that are actually in the container
 ***************************************************************************/
int GTestSuites::test_suites(void) const
{
    return m_testsuites.size();
}

/***********************************************************************//**
 * @brief Return the total number of errors they are in all test suites
 ***************************************************************************/
int GTestSuites::errors(void) const
{
    int errors=0;
    
    //Loop over all test suite
    for(int i=0;i<m_testsuites.size();++i)
    {
        //Add the number of errors
        errors+=m_testsuites[i]->errors(); 
    }
    
    return errors;
}

/***********************************************************************//**
 * @brief Return the total number of failures they are in all test suites
 ***************************************************************************/
int GTestSuites::failures(void) const
{
    int failures=0;
    
    //Loop over all test suite
    for(int i=0;i<m_testsuites.size();++i)
    {
         //Add the number of failures
        failures+=m_testsuites[i]->failures(); 
    }
    
    return failures;
}

/***********************************************************************//**
 * @brief Return the total number of tests they are in all test suites
 ***************************************************************************/
int GTestSuites::tests(void) const
{
    int tests=0;
    
    //Loop over all test suite
    for(int i=0;i<m_testsuites.size();++i)
    {
         //Add the number of tests
        tests+=m_testsuites[i]->tests(); 
    }
    
    return tests; 
}

/***********************************************************************//**
 * @brief Save test repport into XML file.
 *
 * @param[in] filename Name of XML file.
 ***************************************************************************/
void GTestSuites::save(std::string filename)
{
    // Declare empty XML document
    GXml xml;

    // Write observations into XML file
    write(xml);

    // Save XML document
    xml.save(filename);

    // Return
    return;
}

/***********************************************************************//**
 * @brief Return the timestamp. Set at the creation of the object
 ***************************************************************************/
std::time_t GTestSuites::timestamp() const
{
    return m_timestamp;
}

/***********************************************************************//**
 * @brief Write Test Suites into XML document
 *
 * @param[in] xml XML document.
 *
 * Write test suites into the XML document.
 ***************************************************************************/
void GTestSuites::write(GXml& xml)
{
    //Create a xml element for Test Suites
    GXmlElement *element_testsuites = new GXmlElement("testsuites");

    //Set attributes
    element_testsuites->attribute("name",name());
    element_testsuites->attribute("test",str(tests()));
    element_testsuites->attribute("errors",str(errors()));
    element_testsuites->attribute("failures",str(failures()));
    element_testsuites->attribute("time","0"); // not used
    element_testsuites->attribute("timestamp",str(timestamp()));

    //Loop over all test suites in the container
    for (int i=0; i<m_testsuites.size(); ++i)
    {
        //A pointer on the current GTestSuite
        GTestSuite * testsuite = m_testsuites[i];

        //Create a xml element for this test suite
        GXmlElement *element_testsuite = new GXmlElement("testsuite");

        //Set attributes
        element_testsuite->attribute("disabled","");  // not used
        element_testsuite->attribute("errors",str(testsuite->errors()));
        element_testsuite->attribute("failures",str(testsuite->failures()));
        element_testsuite->attribute("hostname","");  // not used
        element_testsuite->attribute("id",str(i));
        element_testsuite->attribute("name",testsuite->name()); 
        element_testsuite->attribute("package","");  // not used
        element_testsuite->attribute("skipped","");  // not used
        element_testsuite->attribute("tests",str(testsuite->tests()));
        element_testsuite->attribute("time","0");  // not used
        element_testsuite->attribute("timestamp",str(testsuite->timestamp()));

        //Loop over all test cases contains in the test suite
        for (int j=0; j<testsuite->tests(); ++j)
        {
            //Reference to the current test case
            GTestCase& testcase=(*testsuite)[j];

            GXmlElement *element_testcase = new GXmlElement("testcase");
            element_testcase->attribute("assertions",""); // not used
            element_testcase->attribute("classname","");  // not used
            element_testcase->attribute("name",testcase.name());
            element_testcase->attribute("status","");  // not used
            element_testcase->attribute("time","0");  // not used

            //If test case is not OK
            if(!testcase.is_passed())
            {
                //Create a xml element for the test case problem.
                GXmlElement* element_testcase_problem = new GXmlElement();

                //Set attributes
                element_testcase_problem->attribute("message",testcase.message());
                element_testcase_problem->attribute("type",testcase.message_type());//TODO not implemented

                //If it is an error
                if(testcase.type()==GTestCase::ERROR_TEST){
                    element_testcase_problem->name("error");
                }
                else{ // else, it is a failure
                    element_testcase_problem->name("failure");
                }
 
                //Append problem to test case element
                element_testcase->append(element_testcase_problem);
            }

            //Append test case to testsuite element.
            element_testsuite->append(element_testcase);
        }

        //Append test suite to test suites element
        element_testsuites->append(element_testsuite);
    }

    //Append test suites to xml document.
    xml.append(element_testsuites);

    //Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                              Private methods                            =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GTestSuites::init_members(){

    //Set default name
    m_name="Unamed Test Suites";

    //Return
    return;
}

/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] testsuites Test suites container.
 *
 * This method just clone the container not the test suite.
 ***************************************************************************/
void GTestSuites::copy_members(const GTestSuites& testsuites){
    m_name = testsuites.m_name;
    m_testsuites = testsuites.m_testsuites;

    return;
}

/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GTestSuites::free_members(){

    return;
}
