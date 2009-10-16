/***************************************************************************
 *                GEbounds.hpp  -  Energy boundary class                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GEbounds.hpp
 * @brief Energy boundary class interface definition.
 * @author J. Knodlseder
 */

#ifndef GBOUNDS_HPP
#define GBOUNDS_HPP

/* __ Includes ___________________________________________________________ */
#include "GFits.hpp"


/***********************************************************************//**
 * @class GEbounds
 *
 * @brief Interface for the GEbounds class.
 ***************************************************************************/
class GEbounds {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GEbounds& ebds);

public:
    // Constructors and destructors
    GEbounds();
    GEbounds(const GEbounds& ebds);
    ~GEbounds();

    // Operators
    GEbounds& operator= (const GEbounds& ebds);

    // Methods
	void        load(const std::string& filename,
                     const std::string& extname = "EBOUNDS");
    void        load(GFitsHDU* hdu);
    int         elements(void) const { return m_num; }
    double      emin(int index) const;
    double      emax(int index) const;
    double      emean(int index) const;
    double      elogmean(int index) const;
    std::string telescope(void) const { return m_telescope; }
    std::string instrument(void) const { return m_instrument; }
    std::string filter(void) const { return m_filter; }
    std::string chantype(void) const { return m_chantype; }
    std::string detchannels(void) const { return m_detchannels; }
  
protected:
    // Protected methods
    void      init_members(void);
    void      copy_members(const GEbounds& ebds);
    void      free_members(void);
    GEbounds* clone(void) const;

    // Protected data area
	int         m_num;           //!< Number of energy boundaries
    int*        m_channel;       //!< Channel number (1,2,3,...)
    double*     m_emin;          //!< Energy bin minima
    double*     m_emax;          //!< Energy bin maxima
    double      m_escale;        //!< Energy scale (1=keV,1000=MeV)
    std::string m_telescope;     //!< Telescope
    std::string m_instrument;    //!< Instrument
    std::string m_filter;        //!< Filter (if any)
    std::string m_chantype;      //!< Channel type
    std::string m_detchannels;   //!< Detector channels
};

#endif /* GBOUNDS_HPP */
