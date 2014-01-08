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
 * This class provides an abstract base class for all FITS table columns.
 * The class supports both fixed-length and variable-length vector columns.
 ***************************************************************************/
class GFitsTableCol : public GBase {

    // Friend classes
    friend class GFitsTable;

public:
    // Constructors and destructors
    GFitsTableCol(void);
    GFitsTableCol(const std::string& name,
                  const int&         length,
                  const int&         number,
                  const int&         width);
    GFitsTableCol(const GFitsTableCol& column);
    virtual ~GFitsTableCol(void);

    // Operators
    GFitsTableCol& operator=(const GFitsTableCol& column);

    // Pure virtual methods
    virtual void            clear(void) = 0;
    virtual GFitsTableCol*  clone(void) const = 0;
    virtual std::string     string(const int& row, const int& inx = 0) const = 0;
    virtual double          real(const int& row, const int& inx = 0) const = 0;
    virtual int             integer(const int& row, const int& inx = 0) const = 0;
    virtual void            insert(const int& row, const int& nrows) = 0;
    virtual void            remove(const int& row, const int& nrows) = 0;
    virtual bool            is_loaded(void) const = 0;

    // Other methods
    void                    name(const std::string& name);
    const std::string&      name(void) const;
    void                    unit(const std::string& unit);
    const std::string&      unit(void) const;
    void                    dim(const std::vector<int>& dim);
    const std::vector<int>& dim(void) const;
    void                    colnum(const int& colnum);
    const int&              colnum(void) const;
    void                    type(const int& type);
    const int&              type(void) const;
    void                    repeat(const int& repeat);
    const int&              repeat(void) const;
    void                    width(const int& width);
    const int&              width(void) const;
    void                    number(const int& number);
    const int&              number(void) const;
    void                    elements(const int& row, const int& elements);
    int                     elements(const int& row) const;
    void                    length(const int& length);
    const int&              length(void) const;
    void                    is_variable(const bool& variable);
    const bool&             is_variable(void) const;
    void                    anynul(const int& anynul);
    const int&              anynul(void) const;
    std::string             tform_binary(void) const;
    std::string             print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void        init_members(void);
    void        copy_members(const GFitsTableCol& column);
    void        free_members(void);
    void        connect(void* vptr);

    // Protected pure virtual methods
    virtual void        alloc_data(void) = 0;
    virtual void        init_data(void) = 0;
    virtual void        fetch_data(void) const = 0;
    virtual void        resize_data(const int& index, const int& number) = 0;
    virtual void        release_data(void) = 0;
    virtual void*       ptr_data(const int& index = 0) = 0;
    virtual void*       ptr_nulval(void) = 0;
    virtual std::string ascii_format(void) const = 0;

    // Protected virtual methods
    virtual void        save(void);
    virtual void        load_column(void);
    virtual void        load_column_fixed(void);
    virtual void        load_column_variable(void);
    virtual void        save_column(void);
    virtual void        save_column_fixed(void);
    virtual void        save_column_variable(void);
    virtual int         offset(const int& row, const int& inx) const;

    // Protected data area
    std::string      m_name;      //!< Column name
    std::string      m_unit;      //!< Column unit
    std::vector<int> m_dim;       //!< Column dimension
    int              m_colnum;    //!< @brief Column number (starting from 1).
                                  //!< This parameter is used to signal if a
                                  //!< table column corresponds to a FITS file
                                  //!< column. If it is set to 0 there is no
                                  //!< correspondance.
    int              m_type;      //!< Column type
    int              m_repeat;    //!< Repeat value of column
    int              m_width;     //!< Width of single column element
    int              m_number;    //!< @brief Number of elements in column.
                                  //!< m_number = m_repeat / m_width
    int              m_length;    //!< Length of column (number of rows)
    bool             m_variable;  //!< Signals if column is variable length
    int              m_varlen;    //!< Maximum number of elements in variable-lgth
    std::vector<int> m_rowstart;  //!< Start index of each row
    mutable int      m_size;      //!< Size of allocated data area (0 if not loaded)
    int              m_anynul;    //!< Number of NULLs encountered
    void*            m_fitsfile;  //!< FITS file pointer associated with column
};


/***********************************************************************//**
 * @brief Set column name
 *
 * @param[in] name Column name.
 ***************************************************************************/
inline
void GFitsTableCol::name(const std::string& name)
{
    // Set name
    m_name = name;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns column name
 *
 * @return Column name.
 ***************************************************************************/
inline
const std::string& GFitsTableCol::name(void) const
{
    // Return name
    return m_name;
}


/***********************************************************************//**
 * @brief Set column unit
 *
 * @param[in] unit Column unit.
 ***************************************************************************/
inline
void GFitsTableCol::unit(const std::string& unit)
{
    // Set unit
    m_unit = unit;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns column unit
 *
 * @return Column unit.
 ***************************************************************************/
inline
const std::string& GFitsTableCol::unit(void) const
{
    // Return column unit
    return m_unit;
}


/***********************************************************************//**
 * @brief Set column dimension
 *
 * @param[in] dim Vector of column dimensions.
 *
 * Sets the column dimension is a integer vector @p dim.
 *
 * @todo Implement dimension check.
 ***************************************************************************/
inline
void GFitsTableCol::dim(const std::vector<int>& dim)
{
    // Set dimension
    m_dim = dim;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns column dimension
 *
 * @return Column dimensions (specified by TDIM keyword).
 ***************************************************************************/
inline
const std::vector<int>& GFitsTableCol::dim(void) const
{
    // Return column dimension
    return m_dim;
}


/***********************************************************************//**
 * @brief Set column number
 *
 * @param[in] colnum Column number.
 ***************************************************************************/
inline
void GFitsTableCol::colnum(const int& colnum)
{
    // Set column number
    m_colnum = colnum;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns column number in FITS file (starting from 1)
 *
 * @return Column number in FITS file (starting from 1).
 ***************************************************************************/
inline
const int& GFitsTableCol::colnum(void) const
{
    // Return column number
    return m_colnum;
}


/***********************************************************************//**
 * @brief Set type code
 *
 * @param[in] type Type code.
 ***************************************************************************/
inline
void GFitsTableCol::type(const int& type)
{
    // Set type code
    m_type = type;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns CFITSIO column type
 *
 * Returns one of the following:
 *   1 (TBIT)
 *  11 (TBYTE)
 *  12 (TSBYTE)
 *  14 (TLOGICAL)
 *  16 (TSTRING)
 *  21 (TUSHORT)
 *  21 (TSHORT)
 *  31 (TUINT)
 *  31 (TINT)
 *  41 (TULONG)
 *  41 (TLONG)
 *  42 (TFLOAT)
 *  81 (TLONGLONG)
 *  82 (TDOUBLE)
 *  83 (TCOMPLEX)
 * 163 (TDBLCOMPLEX)
 *
 * If the type value is negative, the column is a variable-length column.
 ***************************************************************************/
inline
const int& GFitsTableCol::type(void) const
{
    // Return column type
    return m_type;
}


/***********************************************************************//**
 * @brief Set repeat value
 *
 * @param[in] repeat Repeat value.
 ***************************************************************************/
inline
void GFitsTableCol::repeat(const int& repeat)
{
    // Set repeat
    m_repeat = repeat;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns column repeat value (only used for binary tables)
 ***************************************************************************/
inline
const int& GFitsTableCol::repeat(void) const
{
    // Return column repeat value
    return m_repeat;
}


/***********************************************************************//**
 * @brief Set width value
 *
 * @param[in] width Width value.
 ***************************************************************************/
inline
void GFitsTableCol::width(const int& width)
{
    // Set width
    m_width = width;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns width of one column element
 *
 * @return Width of one column element.
 ***************************************************************************/
inline
const int& GFitsTableCol::width(void) const
{
    // Return width of one element in column
    return m_width;
}


/***********************************************************************//**
 * @brief Set number of elements in column
 *
 * @param[in] number Number of elements in column.
 ***************************************************************************/
inline
void GFitsTableCol::number(const int& number)
{
    // Set number of elements
    m_number = number;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns number of elements in column
 *
 * @return Number of elements in column.
 ***************************************************************************/
inline
const int& GFitsTableCol::number(void) const
{
    // Return number of elements in a column
    return m_number;
}


/***********************************************************************//**
 * @brief Set column length (number of rows)
 *
 * @param[in] length Column length.
 *
 * Sets the number of rows in one column.
 ***************************************************************************/
inline
void GFitsTableCol::length(const int& length)
{
    // Set length
    m_length = length;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns number of rows in column
 *
 * @return Number of rows in column.
 ***************************************************************************/
inline
const int& GFitsTableCol::length(void) const
{
    // Return column length
    return m_length;
}


/***********************************************************************//**
 * @brief Set variable-length flag
 *
 * @param[in] variable Variable-length flag.
 ***************************************************************************/
inline
void GFitsTableCol::is_variable(const bool& variable)
{
    // Set variable-length flag
    m_variable = variable;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Signals if column is of variable length
 *
 * @return True if column is a variable length column
 ***************************************************************************/
inline
const bool& GFitsTableCol::is_variable(void) const
{
    // Return variable-length flag
    return m_variable;
}


/***********************************************************************//**
 * @brief Set number of NULLs encountered
 *
 * @param[in] anynul Number of NULLs encountered.
 ***************************************************************************/
inline
void GFitsTableCol::anynul(const int& anynul)
{
    // Set number of NULLs encountered
    m_anynul = anynul;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns number of NULLs encountered
 ***************************************************************************/
inline
const int& GFitsTableCol::anynul(void) const
{
    // Return number of NULLs encountered
    return m_anynul;
}

#endif /* GFITSTABLECOL_HPP */
