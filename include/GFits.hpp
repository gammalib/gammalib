/***************************************************************************
 *                       GFits.hpp - FITS file class                       *
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
 * @file GFits.hpp
 * @brief FITS file class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GFITS_HPP
#define GFITS_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GContainer.hpp"
#include "GFitsHDU.hpp"
#include "GFitsImage.hpp"
#include "GFitsTable.hpp"

/* __ Prototypes _________________________________________________________ */
namespace gammalib {
    int fits_move_to_hdu(const std::string& caller, void* vptr,
                         const int& hdunum = 0);
}


/***********************************************************************//**
 * @class GFits
 *
 * @brief FITS file class
 *
 * This class provides a physical representation of a FITS file in memory.
 * It handles creation, manipulation, writing and reading of FITS files.
 * As a FITS file is a collection of Header Data Units (HDUs), also known
 * as FITS file extensions, the class is designed as a container class of
 * HDUs.
 *
 * HDUs are represented by the abstract GFitsHDU base class. A HDU may be
 * either an image, represented by the abstract GFitsImage base class, or
 * a table, represented by the abstract GFitsTable class. The image() and
 * table() method allow accessing the HDUs. HDUs may be accessed by index
 * (also called extension number) or by extension name.
 ***************************************************************************/
class GFits : public GContainer {

public:
    // Constructors and destructors
    GFits(void);
    explicit GFits(const std::string& filename, const bool& create = false);
    GFits(const GFits& fits);
    virtual ~GFits(void);

    // Operators
    GFits& operator=(const GFits& fits);
    GFitsHDU*       operator[](const int& extno);
    const GFitsHDU* operator[](const int& extno) const;
    GFitsHDU*       operator[](const std::string& extname);
    const GFitsHDU* operator[](const std::string& extname) const;

    // Methods
    void               clear(void);
    GFits*             clone(void) const;
    GFitsHDU*          at(const int& extno);
    const GFitsHDU*    at(const int& extno) const;
    GFitsHDU*          at(const std::string& extname);
    const GFitsHDU*    at(const std::string& extname) const;
    GFitsImage*        image(const int& extno);
    const GFitsImage*  image(const int& extno) const;
    GFitsImage*        image(const std::string& extname);
    const GFitsImage*  image(const std::string& extname) const;
    GFitsTable*        table(const int& extno);
    const GFitsTable*  table(const int& extno) const;
    GFitsTable*        table(const std::string& extname);
    const GFitsTable*  table(const std::string& extname) const;
    int                size(void) const;
    bool               is_empty(void) const;
    GFitsHDU*          set(const int& extno, const GFitsHDU& hdu);
    GFitsHDU*          set(const std::string& extname, const GFitsHDU& hdu);
    GFitsHDU*          append(const GFitsHDU& hdu);
    GFitsHDU*          insert(const int& extno, const GFitsHDU& hdu);
    GFitsHDU*          insert(const std::string& extname, const GFitsHDU& hdu);
    void               remove(const int& extno);
    void               remove(const std::string& extname);
    void               reserve(const int& num);
    void               extend(const GFits& fits);
    bool               contains(const int& extno) const;
    bool               contains(const std::string& extname) const;
    const std::string& filename(void) const;
    int                extno(const std::string& extname) const;
    void               open(const std::string& filename,
                            const bool&        create = false);
    void               save(const bool& clobber = false);
    void               saveto(const std::string& filename,
                              const bool&        clobber = false);
    void               close(void);
    std::string        print(const GChatter& chatter = NORMAL) const;

    // Complex single precision type
    typedef struct {
        float re;
        float im;
    } cfloat;

    // Complex double precision type
    typedef struct {
        float re;
        float im;
    } cdouble;

private:
    // Private methods
    void        init_members(void);
    void        copy_members(const GFits& fits);
    void        free_members(void);
    GFitsImage* new_image(void);
    GFitsImage* new_primary(void);

    // Private data area
    std::vector<GFitsHDU*> m_hdu;        //!< Pointers to HDUs
    std::string            m_filename;   //!< FITS file name
    void*                  m_fitsfile;   //!< FITS file pointer
    bool                   m_readwrite;  //!< FITS file is readwrite (true/false)
    bool                   m_created;    //!< FITS file has been created (true/false)
};


/***********************************************************************//**
 * @brief Get pointer to HDU
 *
 * @param[in] extno Extension number [0,...,size()-1].
 * @return Pointer to HDU.
 *
 * Returns a pointer to the HDU with the specified extension number @p extno.
 * No range checking is performed. If the HDU is not valid, NULL is returned.
 ***************************************************************************/
inline
GFitsHDU* GFits::operator[](const int& extno)
{
    return (m_hdu[extno]);
}


/***********************************************************************//**
 * @brief Get pointer to HDU (const version)
 *
 * @param[in] extno Extension number [0,...,size()-1].
 * @return Pointer to HDU.
 *
 * Returns a pointer to the HDU with the specified extension number @p extno.
 * No range checking is performed. If the HDU is not valid, NULL is returned.
 ***************************************************************************/
inline
const GFitsHDU* GFits::operator[](const int& extno) const
{
    return (m_hdu[extno]);
}


/***********************************************************************//**
 * @brief Get pointer to HDU
 *
 * @param[in] extname Name of HDU extension.
 * @return Pointer to HDU.
 *
 * Returns a pointer to the HDU with the specified @p extname. No checking
 * for the existence of @p extname is performed. If the HDU is not valid,
 * NULL is returned.
 ***************************************************************************/
inline
GFitsHDU* GFits::operator[](const std::string& extname)
{
    return at(extname);
}


/***********************************************************************//**
 * @brief Get pointer to HDU (const version)
 *
 * @param[in] extname Name of HDU extension.
 * @return Pointer to HDU.
 *
 * Returns a pointer to the HDU with the specified @p extname. No checking
 * for the existence of @p extname is performed. If the HDU is not valid,
 * NULL is returned.
 ***************************************************************************/
inline
const GFitsHDU* GFits::operator[](const std::string& extname) const
{
    return at(extname);
}


/***********************************************************************//**
 * @brief Return number of HDUs in FITS file
 *
 * @return Number of HDUs in FITS file.
 *
 * Returns the number of Header Data Units (HDUs) in the FITS file.
 ***************************************************************************/
inline
int GFits::size(void) const
{
    return (m_hdu.size());
}


/***********************************************************************//**
 * @brief Signals if there are no HDUs in FITS file
 *
 * @return True if FITS file is empty, false otherwise.
 *
 * Signals if the FITS file does not contain any HDUs.
 ***************************************************************************/
inline
bool GFits::is_empty(void) const
{
    return (m_hdu.empty());
}


/***********************************************************************//**
 * @brief Reserves space for HDUs in FITS file
 *
 * @param[in] num Number of HDUs
 *
 * Reserves space for @p num HDUs in the FITS file.
 ***************************************************************************/
inline
void GFits::reserve(const int& num)
{
    m_hdu.reserve(num);
    return;
}


/***********************************************************************//**
 * @brief Check if HDU exists in FITS file
 *
 * @param[in] extno Extension number [0,...,size()-1].
 * @return True if HDU with specified @p extno exists, false otherwise.
 *
 * Returns true if a HDU with the specified extension number is present,
 * false otherwise.
 ***************************************************************************/
inline
bool GFits::contains(const int& extno) const
{
    return (extno >= 0 && extno < size());
}


/***********************************************************************//**
 * @brief Check if HDU exists in FITS file
 *
 * @param[in] extname Name of HDU extension.
 * @return True if HDU with specified @p extname exists, false otherwise.
 *
 * Returns true if a HDU with the specified extension name is present,
 * false otherwise.
 ***************************************************************************/
inline
bool GFits::contains(const std::string& extname) const
{
    return (extno(extname) != -1);
}


/***********************************************************************//**
 * @brief Return FITS filename
 *
 * @return FITS filename.
 *
 * Returns the FITS filename. If the object is not yet associated to a file
 * an empty string will be returned.
 ***************************************************************************/
inline
const std::string& GFits::filename(void) const
{
    return (m_filename);
}

#endif /* GFITS_HPP */
