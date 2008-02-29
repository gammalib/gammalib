/***************************************************************************
 *          GSkyDir.hpp  -  Class that implements a sky direction          *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2008 by Jurgen Knodlseder                   *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GSkyDir.hpp
 * @brief GSkyDir class definition.
 * @author J. Knodlseder
 */

#ifndef GSKYDIR_HPP
#define GSKYDIR_HPP

/* __ Includes ___________________________________________________________ */

/* __ Namespaces _________________________________________________________ */


/***************************************************************************
 *                           GSkyDir class definition                      *
 ***************************************************************************/
/**
 * @class GSkyDir
 *
 * @brief Interface for the sky direction classes.
 *
 * @author J. Knodlseder
 */
class GSkyDir {

public:
    /// Constructor
    GSkyDir();

    /// Copy constructor
    GSkyDir(const GSkyDir& dir);

    /// Destructor
    virtual ~GSkyDir();

    /// Assignment operator
    GSkyDir& operator= (const GSkyDir& dir);

protected:
    /// Longitude
    double m_lon;

    /// Latitude
    double m_lat;

    /// Protected method that initialises protected members of the class
    void init_members(void);

    /// Protected method that copies protected members of the class
    void copy_members(const GSkyDir& dir);

    /// Protected method that frees memory allocated by instances of the class
    void free_members(void);

private:

};


/***************************************************************************
 *                              Inline methods                             *
 ***************************************************************************/

#endif /* GSKYDIR_HPP */
