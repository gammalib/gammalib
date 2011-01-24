/***************************************************************************
 *             GEbounds.i  -  Energy boundary class python I/F             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GEbounds.i
 * @brief GEbounds class python interface
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GEbounds.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GEbounds
 *
 * @brief Interface for the GEbounds class.
 *
 * This class holds a list of energy intervals that are used for science
 * analysis.
 ***************************************************************************/
class GEbounds {
public:
    // Constructors and destructors
    GEbounds(void);
    GEbounds(const GEbounds& ebds);
    virtual ~GEbounds(void);

    // Methods
    void    clear(void);
    void    append(const GEnergy& emin, const GEnergy& emax);
    void    insert(const GEnergy& emin, const GEnergy& emax);
    void    setlin(const GEnergy& emin, const GEnergy& emax, const int& num);
    void    setlog(const GEnergy& emin, const GEnergy& emax, const int& num);
	void    load(const std::string& filename,
                 const std::string& extname = "EBOUNDS");
	void    save(const std::string& filename, bool clobber,
                 const std::string& extname = "EBOUNDS") const;
    void    read(GFitsTable* hdu);
    void    write(GFits* file, const std::string& extname = "EBOUNDS") const;
    int     index(const GEnergy& eng) const;
    int     size(void) const { return m_num; }
    GEnergy emin(void) const { return m_emin; }
    GEnergy emax(void) const { return m_emax; }
    GEnergy emin(int inx) const;
    GEnergy emax(int inx) const;
    GEnergy emean(int inx) const;
    GEnergy elogmean(int inx) const;
    bool    isin(const GEnergy& eng) const;
};


/***********************************************************************//**
 * @brief GEbounds class extension
 ***************************************************************************/
%extend GEbounds {
    /*
    char *__str__() {
        return tochar(self->print());
    }
    */
    GEbounds copy() {
        return (*self);
    }
};
