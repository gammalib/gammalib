/***************************************************************************
 *             GHealpix.hpp  -  Healpix sky representation class           *
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
 * @file GHealpix.hpp
 * @brief GHealpix class definition.
 * @author J. Knodlseder
 */

#ifndef GHEALPIX_HPP
#define GHEALPIX_HPP

/* __ Includes ___________________________________________________________ */

/* __ Namespaces _________________________________________________________ */


/***********************************************************************//**
 * @class GHealpix
 *
 * @brief GHealpix class interface defintion
 ***************************************************************************/
class GHealpix {

public:
    // Constructors and destructors
    GHealpix();
    GHealpix(const GHealpix& pixels);
    virtual ~GHealpix();

    // Operators
    GHealpix& operator= (const GHealpix& pixels);

    // Methods
    
protected:
    // Protected methods
    void      init_members(void);
    void      copy_members(const GHealpix& pixels);
    void      free_members(void);
    GHealpix* clone(void) const;

    // Protected data area
private:
};

#endif /* GHEALPIX_HPP */
