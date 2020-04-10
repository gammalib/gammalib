/***************************************************************************
 *     GTestInstDir.hpp  -  Test instrument direction class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2020 by Jean-Baptiste Cayrou                        *
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

#ifndef GTESTINSTDIR_HPP
#define GTESTINSTDIR_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <cstdint>
#include "GInstDir.hpp"


/***********************************************************************//**
 * @class GTestInstDir
 *
 * @brief Test instrument direction class.
 ***************************************************************************/
class GTestInstDir : public GInstDir {

public:
    // Constructors and destructors
    GTestInstDir(void) {
        init_members();
        return;
    }
    GTestInstDir(const GTestInstDir& dir) : GInstDir(dir) {
        init_members();
        copy_members(dir);
        return;
    }
    virtual ~GTestInstDir(void) {
        free_members();
        return;
    }

    // Operators
    GTestInstDir& operator= (const GTestInstDir& dir) {
        if (this != &dir) {
            this->GInstDir::operator=(dir);
            free_members();
            init_members();
            copy_members(dir);
        }
        return *this;
    }

    // Methods
    void clear(void) {
        free_members();
        this->GInstDir::free_members();
        this->GInstDir::init_members();
        init_members();
        return;
    }
    GTestInstDir* clone(void) const {
        return new GTestInstDir(*this);
    }
    std::string classname(void) const {
        return "GTestInstDir";
    }
    void hash(uint64_t hash) {
        m_hash = hash;
    }
    uint64_t hash(void) const {
        return m_hash;
    };
    std::string print(const GChatter& chatter = NORMAL) const {
        std::string result = "=== GTestInstDir ===";
        return result;
    }

protected:
    // Protected methods
    void init_members(void) {
        m_hash = 0;
        return;
    }
    void copy_members(const GTestInstDir& dir) {
        m_hash = dir.m_hash;
        return;
    }
    void free_members(void) {
        return;
    }
    uint64_t m_hash; //!< Hash value
};

#endif /* GTesINSTDIR_HPP */
