/***************************************************************************
 *                GBilinear.hpp - Bilinear interpolator class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2015 by Juergen Knoedlseder                              *
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
 * @file GBilinear.hpp
 * @brief Bilinear interpolator class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GBILINEAR_HPP
#define GBILINEAR_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"


/***********************************************************************//**
 * @class GBilinear
 *
 * @brief Bilinear interpolator class
 ***************************************************************************/
class GBilinear : public GBase {

public:
    // Constructors and destructors
    GBilinear(void);
    GBilinear(const GBilinear& interpolator);
    virtual ~GBilinear(void);

    // Operators
    GBilinear& operator=(const GBilinear& interpolator);
    double     operator()(const double* array);

    // Methods
    void        clear(void);
    GBilinear*  clone(void) const;
    std::string classname(void) const;
    int&        index1(void);
    int&        index2(void);
    int&        index3(void);
    int&        index4(void);
    double&     weight1(void);
    double&     weight2(void);
    double&     weight3(void);
    double&     weight4(void);
    std::string print(const GChatter& chatter = NORMAL) const;

private:
    // Methods
    void init_members(void);
    void copy_members(const GBilinear& interpolator);
    void free_members(void);

    // Members
    int    m_inx1;
    int    m_inx2;
    int    m_inx3;
    int    m_inx4;
    double m_wgt1;
    double m_wgt2;
    double m_wgt3;
    double m_wgt4;
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GBilinear").
 ***************************************************************************/
inline
std::string GBilinear::classname(void) const
{
    return ("GBilinear");
}


/***********************************************************************//**
 * @brief Access index 1
 *
 * @return Reference to index 1.
 ***************************************************************************/
inline
int& GBilinear::index1(void)
{
    return (m_inx1);
}


/***********************************************************************//**
 * @brief Access index 2
 *
 * @return Reference to index 2.
 ***************************************************************************/
inline
int& GBilinear::index2(void)
{
    return (m_inx2);
}


/***********************************************************************//**
 * @brief Access index 3
 *
 * @return Reference to index 3.
 ***************************************************************************/
inline
int& GBilinear::index3(void)
{
    return (m_inx3);
}


/***********************************************************************//**
 * @brief Access index 4
 *
 * @return Reference to index 4.
 ***************************************************************************/
inline
int& GBilinear::index4(void)
{
    return (m_inx4);
}


/***********************************************************************//**
 * @brief Access weight 1
 *
 * @return Reference to weight 1.
 ***************************************************************************/
inline
double& GBilinear::weight1(void)
{
    return (m_wgt1);
}


/***********************************************************************//**
 * @brief Access weight 2
 *
 * @return Reference to weight 2.
 ***************************************************************************/
inline
double& GBilinear::weight2(void)
{
    return (m_wgt2);
}


/***********************************************************************//**
 * @brief Access weight 3
 *
 * @return Reference to weight 3.
 ***************************************************************************/
inline
double& GBilinear::weight3(void)
{
    return (m_wgt3);
}


/***********************************************************************//**
 * @brief Access weight 4
 *
 * @return Reference to weight 4.
 ***************************************************************************/
inline
double& GBilinear::weight4(void)
{
    return (m_wgt4);
}

#endif /* GBILINEAR_HPP */
