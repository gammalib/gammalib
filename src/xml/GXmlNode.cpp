/***************************************************************************
 *             GXmlNode.cpp - XML node base class implementation           *
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
 * @file GXmlNode.cpp
 * @brief XML node base class implementation
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "GException.hpp"
#include "GXmlNode.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_OP_ACCESS1                             "GXmlNode::operator() (int)"
#define G_OP_ACCESS2                       "GXmlNode::operator() (int) const"

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
GXmlNode::GXmlNode(void)
{
    // Initialise private members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] node Object from which the instance should be built.
 ***************************************************************************/
GXmlNode::GXmlNode(const GXmlNode& node)
{
    // Initialise private members for clean destruction
    init_members();

    // Copy members
    copy_members(node);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GXmlNode::~GXmlNode(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Operators                                =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] node Object which should be assigned.
 ***************************************************************************/
GXmlNode& GXmlNode::operator= (const GXmlNode& node)
{
    // Execute only if object is not identical
    if (this != &node) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(node);

    } // endif: object was not identical

    // Return
    return *this;
}


/***********************************************************************//**
 * @brief Node access operator
 *
 * @param[in] index Index of node (0,1,2,...)
 *
 * @exception GException::out_of_range
 *            Node index is out of range.
 ***************************************************************************/
GXmlNode& GXmlNode::operator() (int index)
{
    // If index is outside boundary then throw an error
    if (index < 0 || index >= m_nodes.size())
        throw GException::out_of_range(G_OP_ACCESS1, index, 0, m_nodes.size()-1);

    // Return node
    return *(m_nodes[index]);
}


/***********************************************************************//**
 * @brief Node access operator
 *
 * @param[in] index Index of node (0,1,2,...)
 *
 * @exception GException::out_of_range
 *            Node index is out of range.
 ***************************************************************************/
const GXmlNode& GXmlNode::operator() (int index) const
{
    // If index is outside boundary then throw an error
    if (index < 0 || index >= m_nodes.size())
        throw GException::out_of_range(G_OP_ACCESS2, index, 0, m_nodes.size()-1);

    // Return observation pointer
    return *(m_nodes[index]);
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GXmlNode::init_members(void)
{
    // Initialise members
    m_nodes.clear();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] node Object from which members which should be copied.
 ***************************************************************************/
void GXmlNode::copy_members(const GXmlNode& node)
{
    // Copy nodes
    for (int i = 0; i < node.m_nodes.size(); ++i)
        m_nodes.push_back((node.m_nodes[i]->clone()));

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GXmlNode::free_members(void)
{
    // Free nodes
    for (int i = 0; i < m_nodes.size(); ++i)
        delete m_nodes[i];

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/

