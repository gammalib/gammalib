/***************************************************************************
 *        GFitsTableCol.hpp - FITS table column abstract base class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2013 by Juergen Knoedlseder                         *
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
 * @file GFitsTableCol.hpp
 * @brief FITS table column abstract base class definition
 * @author Juergen Knoedlseder
 */

#ifndef GFITSTABLECOL_HPP
#define GFITSTABLECOL_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GBase.hpp"


/***********************************************************************//**
 * @class GFitsTableCol
 *
 * @brief Abstract interface for FITS table column
 *
 * This class implements a FITS table column. Vector columns are supported.
 ***************************************************************************/
class GFitsTableCol : public GBase {

    // Friend classes
    friend class GFitsTable;

public:
    // Constructors and destructors
    GFitsTableCol(void);
    explicit GFitsTableCol(const std::string& name, const int& length,
                           const int& number, const int& width);
    GFitsTableCol(const GFitsTableCol& column);
    virtual ~GFitsTableCol(void);

    // Operators
    GFitsTableCol& operator= (const GFitsTableCol& column);

    // Pure virtual Methods
    virtual void           clear(void) = 0;
    virtual GFitsTableCol* clone(void) const = 0;
    virtual std::string    string(const int& row, const int& inx = 0) const = 0;
    virtual double         real(const int& row, const int& inx = 0) const = 0;
    virtual int            integer(const int& row, const int& inx = 0) const = 0;
    virtual void           insert(const int& rownum, const int& nrows) = 0;
    virtual void           remove(const int& rownum, const int& nrows) = 0;

    // Base class Methods
    void             name(const std::string& name);
    void             unit(const std::string& unit);
    void             dim(const std::vector<int>& dim);
    std::string      name(void) const;
    std::string      unit(void) const;
    std::vector<int> dim(void) const;
    int              colnum(void) const;
    int              type(void) const;
    int              repeat(void) const;
    int              width(void) const;
    int              number(void) const;
    int              length(void) const;
    int              anynul(void) const;
    std::string      print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected data area
    std::string      m_name;     //!< Column name
    std::string      m_unit;     //!< Column unit
    std::vector<int> m_dim;      //!< Column dimension
    int              m_colnum;   //!< @brief Column number (starting from 1).
                                 //!< This parameter is used to signal if a
                                 //!< table column corresponds to a FITS file
                                 //!< column. If it is set to 0 there is no
                                 //!< correspondance.
    int              m_type;     //!< Column type
    int              m_repeat;   //!< Repeat value of column
    int              m_width;    //!< Width of single column element
    int              m_number;   //!< @brief Number of elements in column.
                                 //!< m_number = m_repeat / m_width
    int              m_length;   //!< Length of column
    mutable int      m_size;     //!< Size of allocated data area (0 if not loaded)
    int              m_anynul;   //!< Number of NULLs encountered
    void*            m_fitsfile; //!< FITS file pointer associated with column

    // Protected pure virtual methods
    virtual std::string ascii_format(void) const = 0;
    virtual std::string binary_format(void) const = 0;
    virtual void        alloc_data(void) = 0;
    virtual void        init_data(void) = 0;
    virtual void*       ptr_data(void) = 0;
    virtual void*       ptr_nulval(void) = 0;

    // Protected virtual methods
    virtual void        save(void);
    virtual void        fetch_data(void) const;
    virtual void        load_column(void);
    virtual void        save_column(void);
    virtual int         offset(const int& row, const int& inx) const;

protected:
    // Private methods
    void init_members(void);
    void copy_members(const GFitsTableCol& column);
    void free_members(void);
    void connect(void* vptr);
};

#endif /* GFITSTABLECOL_HPP */
