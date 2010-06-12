/***************************************************************************
 *                 GGti.hpp  -  Good time interval class                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2010 by Jurgen Knodlseder                           *
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
#include <iostream>
#include "GFits.hpp"
#include "GTime.hpp"


/***********************************************************************//**
 * @class GGti
 *
 * @brief Interface for the GTI class.
 *
 * This class holds a list of Good Time Intervals, i.e. time intervals that
 * are valid for science analysis. Times are stored using the GTime class.
 ***************************************************************************/
class GGti {

	// Friend classes
	friend class GObservation;

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GGti& gti);

public:
    // Constructors and destructors
    GGti(void);
    GGti(const GGti& gti);
    ~GGti(void);

    // Operators
    GGti& operator= (const GGti& gti);

    // Methods
    void   clear(void);
    void   add(const GTime& tstart, const GTime& tstop);
    void   append(const GTime& tstart, const GTime& tstop);
    void   insert(const GTime& tstart, const GTime& tstop);
	void   load(const std::string& filename,
                const std::string& extname = "GTI");
	void   save(const std::string& filename, bool clobber,
                const std::string& extname = "GTI");
    void   read(GFitsHDU* hdu);
    void   write(GFits* file, const std::string& extname = "GTI");
    int    size(void) const { return m_num; }
	GTime  tstart(void) const { return m_tstart; }
	GTime  tstop(void) const { return m_tstop; }
	GTime  tstart(int inx) const;
	GTime  tstop(int inx) const;
	double telapse(void) const { return m_telapse; }
	double ontime(void) const { return m_ontime; }
    bool   isin(const GTime& t) const;
  
protected:
    // Protected methods
    void  init_members(void);
    void  copy_members(const GGti& gti);
    void  free_members(void);
    void  set_attributes(void);
    GGti* clone(void) const;
    void  insert_gti(int inx, const GTime& tstart, const GTime& tstop);
    void  merge_gtis(void);

    // Protected data area
	int     m_num;      //!< Number of intervals
	GTime   m_tstart;   //!< Start of observation
	GTime   m_tstop;    //!< Stop of observation
	double  m_ontime;   //!< Sum of GTI durations (in seconds)
	double  m_telapse;  //!< Time between start of first GTI and stop of last GTI (in seconds)
	GTime  *m_start;    //!< Array of start times
	GTime  *m_stop;     //!< Array of stop times

private:
};

#endif /* GGTI_HPP */
