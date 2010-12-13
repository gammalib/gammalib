/***************************************************************************
 *          GNodeArray.i  -  Array of nodes class SWIG definition          *
 * ----------------------------------------------------------------------- *
 *  copyright : (C) 2008-2010 by Jurgen Knodlseder                         *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GNodeArray.i
 * @brief GNodeArray class SWIG file.
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GNodeArray.hpp"
%}


/***********************************************************************//**
 * @class GNodeArray
 *
 * @brief SWIG interface for the node array class.
 ***************************************************************************/
class GNodeArray {
public:
    // Constructors and destructors
    GNodeArray(void);
    GNodeArray(const GNodeArray& array);
    virtual ~GNodeArray(void);

    // Methods
    void        clear(void);
    GNodeArray* clone(void) const;
    int         size(void) { return m_node.size(); }
    void        nodes(const int& num, const double* array);
    void        nodes(const GVector& vector);
    void        nodes(const std::vector<double>& vector);
    void        append(const double& node);
    double      interpolate(const double& value, const std::vector<double>& vector);
    void        set_value(const double& value);
    int         inx_left(void) { return m_inx_left; }
    int         inx_right(void) { return m_inx_right; }
    double      wgt_left(void) { return m_wgt_left; }
    double      wgt_right(void) { return m_wgt_right; }
};


/***********************************************************************//**
 * @brief GNodeArray class extension
 ***************************************************************************/
%extend GNodeArray {
    /*
    char *__str__() {
        static std::string result = self->print();
        return ((char*)result.c_str());
    }
    */    
    GNodeArray copy() {
        return (*self);
    }
};
