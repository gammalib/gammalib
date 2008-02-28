/***************************************************************************
 *        GFitsTableCol.hpp  - FITS table column abstract base class       *
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

#ifndef GFITSTABLECOL_HPP
#define GFITSTABLECOL_HPP

/* __ Includes ___________________________________________________________ */
#include "GFitsCfitsio.hpp"

/* __ Namespaces _________________________________________________________ */


/***************************************************************************
 *                      GFitsTableCol class definition                     *
 ***************************************************************************/
class GFitsTableCol {

public:
    // Constructors and destructors
    GFitsTableCol();
    GFitsTableCol(const GFitsTableCol& column);
    virtual ~GFitsTableCol();

    // Operators
    GFitsTableCol& operator= (const GFitsTableCol& column);

    // Virtual Methods
    virtual std::string    string(const int& row, const int& col = 0) = 0;
    virtual double         real(const int& row, const int& col = 0) = 0;
    virtual int            integer(const int& row, const int& col = 0) = 0;
    virtual GFitsTableCol* clone(void) const = 0;
    virtual std::string*   ptr_string(void) = 0;
    virtual double*        ptr_double(void) = 0;
    virtual float*         ptr_float(void) = 0;
    virtual short*         ptr_short(void) = 0;
    virtual long*          ptr_long(void) = 0;
    virtual int*           ptr_int(void) = 0;

    // Base class Methods
    void        set_name(const std::string name);
    void        set_colnum(const int colnum);
    void        set_type(const int type);
    void        set_repeat(const int repeat);
    void        set_width(const int width);
    void        set_length(const int length);
    void        set_fitsfile(const __fitsfile* fptr);
    std::string name(void);
    int         colnum(void);
    int         type(void);
    int         repeat(void);
    int         width(void);
    int         length(void);
    
protected:
    // Protected data area
    std::string m_name;
    int         m_colnum;
    int         m_type;
    int         m_repeat;
    int         m_width;
    int         m_length;
    __fitsfile* m_fitsfile;

private:
    // Private methods
    void init_members(void);
    void copy_members(const GFitsTableCol& column);
    void free_members(void);
};


/***************************************************************************
 *                              Inline methods                             *
 ***************************************************************************/
inline void         GFitsTableCol::set_name(const std::string name) { m_name = name; }
inline void         GFitsTableCol::set_colnum(const int colnum) { m_colnum = colnum; }
inline void         GFitsTableCol::set_type(const int type) { m_type = type; }
inline void         GFitsTableCol::set_repeat(const int repeat) { m_repeat = repeat; }
inline void         GFitsTableCol::set_width(const int width) { m_width = width; }
inline void         GFitsTableCol::set_length(const int length) { m_length = length; }
inline void         GFitsTableCol::set_fitsfile(const __fitsfile* fptr) { m_fitsfile = (__fitsfile*)fptr; }
inline std::string  GFitsTableCol::name(void) { return m_name; }
inline int          GFitsTableCol::colnum(void) { return m_colnum; }
inline int          GFitsTableCol::type(void) { return m_type; }
inline int          GFitsTableCol::repeat(void) { return m_repeat; }
inline int          GFitsTableCol::width(void) { return m_width; }
inline int          GFitsTableCol::length(void) { return m_length; }

#endif /* GFITSTABLECOL_HPP */
