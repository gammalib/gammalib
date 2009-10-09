/***************************************************************************
 *                 GGti.hpp  -  Good time interval class                   *
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
 * @file GGti.hpp
 * @brief Good time interval class interface definition.
 * @author J. Knodlseder
 */

#ifndef GGTI_HPP
#define GGTI_HPP

/* __ Includes ___________________________________________________________ */


/***********************************************************************//**
 * @class GGti
 *
 * @brief Interface for the GTI class.
 ***************************************************************************/
class GGti {

  // Friend classes
  friend class GObservation;

public:
    // Constructors and destructors
    GGti();
    GGti(const GGti& gti);
    ~GGti();

    // Operators
    GGti& operator= (const GGti& gti);

    // Methods
  
protected:
    // Protected methods
    void  init_members(void);
    void  copy_members(const GGti& gti);
    void  free_members(void);
    GGti* clone(void) const;

    // Protected data area
	int         m_num;          //!< Number of intervals
	double      m_tstart;       //!< Start of observation
	double      m_tstop;        //!< Stop of observation
	double      m_ontime;       //!< Sum of GTI durations
	double      m_elapse;       //!< Time between start of first GTI and stop of last GTI
	double     *m_start;        //!< Array of start times
	double     *m_stop;         //!< Array of stop times

private:
};

#endif /* GGTI_HPP */
