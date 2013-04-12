/***************************************************************************
 *       GFitsHeaderCard.cpp  - FITS header card abstract base class       *
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
 * @file GFitsHeaderCard.cpp
 * @brief FITS header card class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cfloat>
#include <climits>
#include "GException.hpp"
#include "GTools.hpp"
#include "GFitsCfitsio.hpp"
#include "GFitsHeaderCard.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_READ_NUM                         "GFitsHeaderCard::read(void*,int)"
#define G_READ_STR                 "GFitsHeaderCard::read(void*,std::string)"
#define G_WRITE                               "GFitsHeaderCard::write(void*)"

/* __ Definitions ________________________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Enumerations _______________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GFitsHeaderCard::GFitsHeaderCard(void)
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor for string cards
 *
 * @param[in] keyname Card name.
 * @param[in] value Card string value.
 * @param[in] comment Card comment.
 *
 * This constructor builds a header card from the keyname, value and comment.
 ***************************************************************************/
GFitsHeaderCard::GFitsHeaderCard(const std::string& keyname,
                                 const std::string& value,
                                 const std::string& comment)
{
    // Initialise class members for clean destruction
    init_members();

    // Set members
    this->keyname(keyname);
    this->value("'"+value+"'");
    this->comment(comment);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor for floating point cards
 *
 * @param[in] keyname Card name.
 * @param[in] value Card floating point value.
 * @param[in] comment Card comment.
 *
 * This constructor builds a header card from the keyname, value and comment.
 ***************************************************************************/
GFitsHeaderCard::GFitsHeaderCard(const std::string& keyname,
                                 const double&      value,
                                 const std::string& comment)
{
    // Initialise class members for clean destruction
    init_members();

    // Set members
    this->keyname(keyname);
    this->value(value);
    this->comment(comment);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor for integer cards
 *
 * @param[in] keyname Card name.
 * @param[in] value Card integer value.
 * @param[in] comment Card comment.
 *
 * This constructor builds a header card from the keyname, value and comment.
 ***************************************************************************/
GFitsHeaderCard::GFitsHeaderCard(const std::string& keyname,
                                 const int&         value,
                                 const std::string& comment)
{
    // Initialise class members for clean destruction
    init_members();

    // Set members
    this->keyname(keyname);
    this->value(value);
    this->comment(comment);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] card Header card.
 ***************************************************************************/
GFitsHeaderCard::GFitsHeaderCard(const GFitsHeaderCard& card)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(card);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GFitsHeaderCard::~GFitsHeaderCard(void)
{
    // Free members
    free_members();

    // Return
    return;
}

/*==========================================================================
 =                                                                         =
 =                                Operators                                =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] card Header card.
 ***************************************************************************/
GFitsHeaderCard& GFitsHeaderCard::operator= (const GFitsHeaderCard& card)
{
    // Execute only if object is not identical
    if (this != &card) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(card);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear object
 ***************************************************************************/
void GFitsHeaderCard::clear(void)
{
    // Free members
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone object
 ***************************************************************************/
GFitsHeaderCard* GFitsHeaderCard::clone(void) const
{
    return new GFitsHeaderCard(*this);
}


/***********************************************************************//**
 * @brief Set name of header card
 *
 * @param[in] keyname Card name.
 ***************************************************************************/
void GFitsHeaderCard::keyname(const std::string& keyname)
{
    // Set name of card
    m_keyname = keyname;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set string value of header card
 *
 * @param[in] value Card string value.
 *
 * This method sets the value of a string header card. The internal string
 * representation contains hyphens, yet for keyword writing the hyphens are
 * stripped.
 ***************************************************************************/
void GFitsHeaderCard::value(const std::string& value)
{
    // Free data type
    free_dtype();

    // Set value
    m_value = value;

    // Attach hyphens to internal string representation if required
    if (m_value[0] != '\x27')
        m_value = "'" + m_value;
    if (m_value[m_value.length()-1] != '\x27')
        m_value = m_value + "'";

    // Strip hyphens and whitespace from datatype value that is used for
    // keyword writing
    std::string value_dtype = 
                gammalib::strip_whitespace(m_value.substr(1, m_value.length() - 2));

    // Set data type
    m_dtype       = __TSTRING;
    m_value_dtype = new std::string(value_dtype);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set boolean value of header card
 *
 * @param[in] value Card boolean value.
 ***************************************************************************/
void GFitsHeaderCard::value(const bool& value)
{
    // Free data type
    free_dtype();

    // Set value and data type
    m_value       = (value) ? "T" : "F";
    m_dtype       = __TLOGICAL;
    m_value_dtype = new bool(value);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set single precision value of header card
 *
 * @param[in] value Card single precision value.
 ***************************************************************************/
void GFitsHeaderCard::value(const float& value)
{
    // Free data type
    free_dtype();

    // Set value and data type
    m_value       = gammalib::str(value);
    m_dtype       = __TFLOAT;
    m_value_dtype = new float(value);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set double precision value of header card
 *
 * @param[in] value Card double precision value.
 ***************************************************************************/
void GFitsHeaderCard::value(const double& value)
{
    // Free data type
    free_dtype();

    // Set value and data type
    m_value       = gammalib::str(value);
    m_dtype       = __TDOUBLE;
    m_value_dtype = new double(value);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set unsigned short integer value of header card
 *
 * @param[in] value Card unsigned short integer value.
 ***************************************************************************/
void GFitsHeaderCard::value(const unsigned short& value)
{
    // Free data type
    free_dtype();

    // Set value and data type
    m_value       = gammalib::str(value);
    m_dtype       = __TUSHORT;
    m_value_dtype = new unsigned short(value);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set short integer value of header card
 *
 * @param[in] value Card short integer value.
 ***************************************************************************/
void GFitsHeaderCard::value(const short& value)
{
    // Free data type
    free_dtype();

    // Set value and data type
    m_value       = gammalib::str(value);
    m_dtype       = __TSHORT;
    m_value_dtype = new short(value);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set unsigned integer value of header card
 *
 * @param[in] value Card unsigned integer value.
 ***************************************************************************/
void GFitsHeaderCard::value(const unsigned int& value)
{
    // Free data type
    free_dtype();

    // Set value and data type
    m_value       = gammalib::str(value);
    m_dtype       = __TUINT;
    m_value_dtype = new unsigned int(value);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set integer value of header card
 *
 * @param[in] value Card integer value.
 ***************************************************************************/
void GFitsHeaderCard::value(const int& value)
{
    // Free data type
    free_dtype();

    // Set value and data type
    m_value       = gammalib::str(value);
    m_dtype       = __TINT;
    m_value_dtype = new int(value);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set long integer value of header card
 *
 * @param[in] value Card long integer value.
 ***************************************************************************/
void GFitsHeaderCard::value(const long& value)
{
    // Free data type
    free_dtype();

    // Set value and data type
    m_value       = gammalib::str(value);
    m_dtype       = __TLONG;
    m_value_dtype = new long(value);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set unsigned integer long value of header card
 *
 * @param[in] value Card unsigned long integer value.
 ***************************************************************************/
void GFitsHeaderCard::value(const unsigned long& value)
{
    // Free data type
    free_dtype();

    // Set value and data type
    m_value       = gammalib::str(value);
    m_dtype       = __TULONG;
    m_value_dtype = new unsigned long(value);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set long long integer value of header card
 *
 * @param[in] value Card long long integer value.
 ***************************************************************************/
void GFitsHeaderCard::value(const long long& value)
{
    // Free data type
    free_dtype();

    // Set value and data type
    m_value       = gammalib::str(value);
    m_dtype       = __TLONGLONG;
    m_value_dtype = new long long(value);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set unit of header card value
 *
 * @param[in] unit Unit of header card.
 ***************************************************************************/
void GFitsHeaderCard::unit(const std::string& unit)
{
    // Set unit
    m_unit = unit;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set comment of header card
 *
 * @param[in] comment Header card comment.
 ***************************************************************************/
void GFitsHeaderCard::comment(const std::string& comment)
{
    // Set comment
    m_comment = comment;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Return header card keyname
  ***************************************************************************/
std::string GFitsHeaderCard::keyname(void) const
{
    // Return
    return m_keyname;
}


/***********************************************************************//**
 * @brief Return header card value
  ***************************************************************************/
std::string GFitsHeaderCard::value(void) const
{
    // Return
    return m_value;
}


/***********************************************************************//**
 * @brief Return header card decimals
  ***************************************************************************/
int GFitsHeaderCard::decimals(void) const
{
    // Return
    return m_value_decimals;
}


/***********************************************************************//**
 * @brief Return header card value unit
  ***************************************************************************/
std::string GFitsHeaderCard::unit(void) const
{
    // Return
    return m_unit;
}


/***********************************************************************//**
 * @brief Return header card comment
  ***************************************************************************/
std::string GFitsHeaderCard::comment(void) const
{
    // Return
    return m_comment;
}


/***********************************************************************//**
 * @brief Return header card value as string
 *
 * Convert header card value into a string.
 * Any hyphens that may occur in the FITS card will be automatically stripped.
 ***************************************************************************/
std::string GFitsHeaderCard::string(void)
{
    // Initialize return value to actual value string
    std::string result = m_value;

    // Type dependent conversion
    if (m_value_dtype != NULL) {
        switch (m_dtype) {
        case __TSTRING:
            if (m_value.length() > 2)
                result = gammalib::strip_whitespace(m_value.substr(1, m_value.length() - 2));
            break;
        default:
            break;
        }
    }

    // Return string
    return result;

}


/***********************************************************************//**
 * @brief Return header card value as double precision
 *
 * Convert header card value into a double precision value. In case that the
 * card did not contain a numerical value, 0 will be returned by the
 * method.
 ***************************************************************************/
double GFitsHeaderCard::real(void)
{
    // Initialize return value to 0.0
    double result = 0.0;

    // Type dependent conversion
    if (m_value_dtype != NULL) {
        switch (m_dtype) {
        case __TSTRING:
            if (m_value.length() > 2)
                result = gammalib::todouble(m_value.substr(1, m_value.length() - 2));
            break;
        case __TLOGICAL:
            result = (m_value == "T") ? 1.0 : 0.0;
            break;
        default:
            result = gammalib::todouble(m_value);
            break;
        }
    }

    // Return string
    return result;

}


/***********************************************************************//**
 * @brief Return header card value as integer
*
 * Convert header card value into a integer value. In case that the
 * card did not contain a numerical value, 0 will be returned by the
 * method.
 ***************************************************************************/
int GFitsHeaderCard::integer(void)
{
    // Initialize return value to 0
    int result = 0;

    // Type dependent conversion
    if (m_value_dtype != NULL) {
        switch (m_dtype) {
        case __TSTRING:
            if (m_value.length() > 2)
                result = gammalib::toint(m_value.substr(1, m_value.length() - 2));
            break;
        case __TLOGICAL:
            result = (m_value == "T") ? 1 : 0;
            break;
        default:
            result = gammalib::toint(m_value);
            break;
        }
    }

    // Return string
    return result;

}


/***********************************************************************//**
 * @brief Print header card information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing header card information.
 ***************************************************************************/
std::string GFitsHeaderCard::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;
    
    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append keyname
        result.append(gammalib::left(m_keyname,8));

        // Format values
        if (m_keyname != "COMMENT" && m_keyname != "HISTORY") {
            if (m_unit.length() > 0) {
                result.append(" ="+gammalib::right(m_value,21)+" / ["+m_unit+"] "+m_comment);
            }
            else {
                result.append(" ="+gammalib::right(m_value,21)+" / "+m_comment);
            }
        }
        else {
            result.append(" "+m_comment);
        }

    // Attach card type
    if (m_value_dtype != NULL) {
        switch (m_dtype) {
        case __TBIT:
            result.append(" <bit>");
            break;
        case __TBYTE:
            result.append(" <byte>");
            break;
        case __TSBYTE:
            result.append(" <signed byte>");
            break;
        case __TLOGICAL:
            result.append(" <bool>");
            break;
        case __TSTRING:
            result.append(" <string>");
            break;
        case __TUSHORT:
            result.append(" <unsigned short>");
            break;
        case __TSHORT:
            result.append(" <short>");
            break;
        case __TUINT:
            result.append(" <unsigned int>");
            break;
        case __TINT:
            result.append(" <int>");
            break;
        case __TULONG:
            result.append(" <unsigned long>");
            break;
        case __TLONG:
            result.append(" <long>");
            break;
        case __TLONGLONG:
            result.append(" <long long>");
            break;
        case __TFLOAT:
            result.append(" <float>");
            break;
        case __TDOUBLE:
            result.append(" <double>");
            break;
        case __TCOMPLEX:
            result.append(" <complex>");
            break;
        case __TDBLCOMPLEX:
            result.append(" <double complex>");
            break;
        default:
            result.append(" <unsupported>");
            break;
        }
    }
    else {
        result.append(" <non native>");
    }

    } // endif: chatter was not silent

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GFitsHeaderCard::init_members(void)
{
    // Initialise members
    m_keyname.clear();
    m_value.clear();
    m_unit.clear();
    m_comment.clear();
    m_value_dtype    = NULL;
    m_dtype          = 0;
    m_value_decimals = 10;
    m_comment_write  = false;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] card Header card to be copied
 ***************************************************************************/
void GFitsHeaderCard::copy_members(const GFitsHeaderCard& card)
{
    // Copy members
    m_keyname        = card.m_keyname;
    m_value          = card.m_value;
    m_value_decimals = card.m_value_decimals;
    m_unit           = card.m_unit;
    m_comment        = card.m_comment;
    m_comment_write  = card.m_comment_write;

    // Copy native data types
    copy_dtype(card);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GFitsHeaderCard::free_members(void)
{
    // Free members
    free_dtype();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy dtype
 ***************************************************************************/
void GFitsHeaderCard::copy_dtype(const GFitsHeaderCard& card)
{
    // Copy data type
    if (card.m_value_dtype != NULL) {
        m_dtype = card.m_dtype;
        switch (m_dtype) {
        case __TLOGICAL:
            m_value_dtype = new bool(*((bool*)card.m_value_dtype));
            break;
        case __TSTRING:
            m_value_dtype = new std::string(*((std::string*)card.m_value_dtype));
            break;
        case __TUSHORT:
            m_value_dtype = new unsigned short(*((unsigned short*)card.m_value_dtype));
            break;
        case __TSHORT:
            m_value_dtype = new short(*((short*)card.m_value_dtype));
            break;
        case __TUINT:
            m_value_dtype = new unsigned int(*((unsigned int*)card.m_value_dtype));
            break;
        case __TINT:
            m_value_dtype = new int(*((int*)card.m_value_dtype));
            break;
        case __TULONG:
            m_value_dtype = new unsigned long(*((unsigned long*)card.m_value_dtype));
            break;
        case __TLONG:
            m_value_dtype = new long(*((long*)card.m_value_dtype));
            break;
        case __TLONGLONG:
            m_value_dtype = new long long(*((long long*)card.m_value_dtype));
            break;
        case __TFLOAT:
            m_value_dtype = new float(*((float*)card.m_value_dtype));
            break;
        case __TDOUBLE:
            m_value_dtype = new double(*((double*)card.m_value_dtype));
            break;
        default:
            std::cout << "GFitsHeaderCard::copy_dtype: invalid data type code "
                      << m_dtype << " encountered." << std::endl;
            break;
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Free dtype
 ***************************************************************************/
void GFitsHeaderCard::free_dtype(void)
{
    // Free data type
    if (m_value_dtype != NULL) {
        switch (m_dtype) {
        case __TLOGICAL:
            delete (bool*)m_value_dtype;
            break;
        case __TSTRING:
            delete (std::string*)m_value_dtype;
            break;
        case __TUSHORT:
            delete (unsigned short*)m_value_dtype;
            break;
        case __TSHORT:
            delete (short*)m_value_dtype;
            break;
        case __TUINT:
            delete (unsigned int*)m_value_dtype;
            break;
        case __TINT:
            delete (int*)m_value_dtype;
            break;
        case __TULONG:
            delete (unsigned long*)m_value_dtype;
            break;
        case __TLONG:
            delete (long*)m_value_dtype;
            break;
        case __TLONGLONG:
            delete (long long*)m_value_dtype;
            break;
        case __TFLOAT:
            delete (float*)m_value_dtype;
            break;
        case __TDOUBLE:
            delete (double*)m_value_dtype;
            break;
        default:
            std::cout << "GFitsHeaderCard::free_dtype: invalid data type code "
                      << m_dtype << " encountered." << std::endl;
            break;
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set native data type from card string
 *
 * @param[in] value Card string.
 *
 * The native data type is determined from the card string. This method sets
 * the members m_dtype and m_value_type. The card type is determined by
 * analyzing the card value.
 *
 * @todo Implement syntax checking for integer or float values.
 ***************************************************************************/
void GFitsHeaderCard::set_dtype(const std::string& value)
{
    // Free data type
    free_dtype();

    // Main loop to allow fall through
    do {

        // Get index of last string element
        int last = m_value.length() - 1;

        // If value is empty then we either have a COMMENT or HISTORY card or
        // we don't know the type
        if (last < 0) {
            if (m_keyname == "COMMENT") {
                continue;
            }
            else if (m_keyname == "HISTORY") {
                continue;
            }
            else {
                // TODO: If we reach this point we have an unrecognised keyname
                continue;
            }
        }

        // If value is enclosed in parentheses then we have a string
        if (m_value[0] == '\x27' && m_value[last] == '\x27') {
            this->value(value);
            continue;
        }

        // If value has only one digit and is either 'F' or 'T' we have a
        // Boolean
        if (last == 0 && (m_value[0] == 'F' || m_value[0] == 'T')) {
            this->value((m_value[0] == 'T'));
            continue;
        }

        // Conversion
        double             value_dbl = gammalib::todouble(m_value);
        long long          value_ll  = gammalib::tolonglong(m_value);
        unsigned long long value_ull = gammalib::toulonglong(m_value);

        // Check if we have an integer
        if ((value_dbl >= 0 && value_dbl == value_ull) ||
            (value_dbl <  0 && value_dbl == value_ll)) {

            // Consider positive integers as unsigned
            if (value_dbl >= 0) {
                if (value_ull > ULONG_MAX) {
                    m_dtype       = __TLONGLONG;
                    m_value_dtype = new long long(value_ll);
                }
                else if (value_ull > USHRT_MAX) {
                    m_dtype       = __TULONG;
                    m_value_dtype = new unsigned long(value_ull);
                }
                else {
                    m_dtype       = __TUSHORT;
                    m_value_dtype = new unsigned short(value_ull);
                }
            }

            // ... otherwise consider signed
            else {
                if (value_ll > LONG_MAX || value_ll < LONG_MIN) {
                    m_dtype       = __TLONGLONG;
                    m_value_dtype = new long long(value_ll);
                }
                else if (value_ll > SHRT_MAX || value_ll < SHRT_MIN) {
                    m_dtype       = __TLONG;
                    m_value_dtype = new long(value_ll);
                }
                else {
                    m_dtype       = __TSHORT;
                    m_value_dtype = new short(value_ll);
                }
            }
        } // endif: handled integers

        // ... otherwise handle floats
        else {
            if (value_dbl > FLT_MAX || value_dbl < FLT_MIN) {
                m_dtype       = __TDOUBLE;
                m_value_dtype = new double(value_dbl);
            }
            else {
                m_dtype       = __TFLOAT;
                m_value_dtype = new float(value_dbl);
            }
        }

    } while(0);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read header card from FITS file
 *
 * @param[in] vptr FITS file void pointer.
 * @param[in] keynum Number of the header card.
 ***************************************************************************/
void GFitsHeaderCard::read(void* vptr, int keynum)
{
    // Move to HDU
    int status = 0;
    status     = __ffmahd(FPTR(vptr), (FPTR(vptr)->HDUposition)+1, NULL,
                          &status);
    if (status != 0)
        throw GException::fits_error(G_READ_NUM, status);

    // Read keyword
    char keyname[80];
    char value[80];
    char comment[80];
    status = __ffgkyn(FPTR(vptr), keynum, keyname, value, comment, &status);
    if (status != 0)
        throw GException::fits_error(G_READ_NUM, status);

    // Store result
    m_keyname.assign(keyname);
    m_value.assign(value);
    m_comment.assign(comment);

    // Determine card type
    set_dtype(m_value);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Read header card from FITS file
 *
 * @param[in] vptr FITS file void pointer.
 * @param[in] keyname Name of the header card.
 ***************************************************************************/
void GFitsHeaderCard::read(void* vptr, const std::string& keyname)
{
    // Move to HDU
    int status = 0;
    status     = __ffmahd(FPTR(vptr), (FPTR(vptr)->HDUposition)+1, NULL,
                          &status);
    if (status != 0)
        throw GException::fits_error(G_READ_NUM, status);

    // Read keyword
    char value[80];
    char comment[80];
    status = __ffgkey(FPTR(vptr), (char*)keyname.c_str(), value, comment,
                      &status);

    // Catch error
    if (status == 202)       // Keyword not found
        throw GException::fits_key_not_found(G_READ_STR, keyname, status);
    else if (status != 0)    // Any other error
        throw GException::fits_error(G_READ_STR, status);

    // Store result
    m_keyname = keyname;
    m_value.assign(value);
    m_comment.assign(comment);

    // Determine card type
    set_dtype(m_value);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Write header card
 *
 * @param[in] vptr FITS file void pointer.
 *
 * Writes any kind of header card to a FITS file.
 ***************************************************************************/
void GFitsHeaderCard::write(void* vptr)
{
    // Move to HDU
    int status = 0;
    status     = __ffmahd(FPTR(vptr), (FPTR(vptr)->HDUposition)+1, NULL,
                          &status);
    if (status != 0)
        throw GException::fits_error(G_READ_NUM, status);

    // If card is comment then write comment
    if (m_keyname == "COMMENT") {
        if (m_comment_write)
            status = __ffpcom(FPTR(vptr), (char*)m_comment.c_str(), &status);
    }

    // If card is history then write history
    else if (m_keyname == "HISTORY") {
        if (m_comment_write)
            status = __ffphis(FPTR(vptr), (char*)m_comment.c_str(), &status);
    }

    // If card holds a native data type then write it
    else if (m_value_dtype != NULL) {
        switch (m_dtype) {
        case __TLOGICAL:
            {
            int value = (*((bool*)m_value_dtype)) ? 1 : 0;
            status = __ffukyl(FPTR(vptr), (char*)m_keyname.c_str(), value,
                              (char*)m_comment.c_str(), &status);
            }
            break;
        case __TSTRING:
            status = __ffukys(FPTR(vptr), (char*)m_keyname.c_str(),
                              (char*)((std::string*)m_value_dtype)->c_str(),
                              (char*)m_comment.c_str(), &status);
            break;
        case __TFLOAT:
            status = __ffukye(FPTR(vptr), (char*)m_keyname.c_str(),
                              *((float*)m_value_dtype), decimals(),
                              (char*)m_comment.c_str(), &status);
            break;
        case __TDOUBLE:
            status = __ffukyd(FPTR(vptr), (char*)m_keyname.c_str(),
                              *((double*)m_value_dtype), decimals(),
                              (char*)m_comment.c_str(), &status);
            break;
        default:
            status = __ffuky(FPTR(vptr), m_dtype, (char*)m_keyname.c_str(),
                             m_value_dtype, (char*)m_comment.c_str(), &status);
            break;
        }
    } // endif: had native data type

    // ... capture all other stuff
    else {
        std::cout << "*** ERROR in GFitsHeaderCard::write: the code should"
                  << " never arrive at this point. keyname=" << m_keyname
                  << " " << "dtype=" << m_dtype << std::endl;
    }

    // Throw exception in case of a FITS error
    if (status != 0)
        throw GException::fits_error(G_WRITE, status);

    // Return
    return;
}
