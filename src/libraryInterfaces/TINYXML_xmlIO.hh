/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2015 FlowKit Sarl
 * Route d'Oron 2
 * 1010 Lausanne, Switzerland
 * E-mail contact: contact@flowkit.com
 *
 * The most recent release of Palabos can be downloaded at 
 * <http://www.palabos.org/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/** \file
 * Input/Output in XML format -- generic code.
 */

#ifndef XML_IO_HH
#define XML_IO_HH

#include "libraryInterfaces/TINYXML_xmlIO.h"
#include <typeinfo>
#include <cctype>
#include <algorithm>
#include <string> 

namespace epc {

template <typename T>
void XMLreaderProxy::read(T& value) const {
    if (!reader) return;
    std::stringstream valueStr(reader->getText(id));
    T tmp = T();
    if (!(valueStr>>tmp)) {
        EPC_ASSERT(false && (std::string("Cannot read value from XML element ") + reader->getName()).c_str());
    }
    value = tmp;
}

template <>
inline void XMLreaderProxy::read<bool>(bool& value) const {
    if (!reader) return;
    std::stringstream valueStr(reader->getText(id));
    std::string word;
    valueStr >> word;
    // Transform to lower-case, so that "true" and "false" are case-insensitive.
    std::transform(word.begin(), word.end(), word.begin(), ::tolower);
    if (word=="true") {
        value = true;
    }
    else if (word=="false") {
        value=false;
    }
    else {
        EPC_ASSERT(false && (std::string("Cannot read value from XML element ") + reader->getName()).c_str());
    }
}

template <>
inline void XMLreaderProxy::read<std::string>(std::string& entry) const {
    if (!reader) return;
    entry = reader->getText(id);
}

template <typename T>
bool XMLreaderProxy::readNoThrow(T& value) const {
    if (!reader) return false;
    std::stringstream valueStr(reader->getText(id));
    T tmp = T();
    if (!(valueStr >> tmp)) {
        return false;
    }
    value = tmp;
    return true;
}

template <>
inline bool XMLreaderProxy::readNoThrow<bool>(bool& value) const {
    if (!reader) return false;
    std::stringstream valueStr(reader->getText(id));
    std::string word;
    valueStr >> word;
    // Transform to lower-case, so that "true" and "false" are case-insensitive.
    std::transform(word.begin(), word.end(), word.begin(), ::tolower);
    if (word=="true") {
        value = true;
        return true;
    }
    else if (word=="false") {
        value=false;
        return true;
    }
    return false;
}

template <>
inline bool XMLreaderProxy::readNoThrow<std::string>(std::string& entry) const {
    if (!reader) return false;
    entry = reader->getText(id);
    return true;
}

template <typename T>
void XMLreaderProxy::read(std::vector<T>& values) const {
    if (!reader) return;
    std::stringstream multiValueStr(reader->getText(id));
    std::string word;
    std::vector<T> tmp(values);
    while (multiValueStr>>word) {
        std::stringstream valueStr(word);
        T value;
        if (!(valueStr >> value)) {
            EPC_ASSERT(false && (std::string("Cannot read value from XML element ") + reader->getName()).c_str());
        }
        tmp.push_back(value);
    }
    values.swap(tmp);
}

template <typename T>
bool XMLreaderProxy::readNoThrow(std::vector<T>& values) const {
    if (!reader) return false;
    std::stringstream multiValueStr(reader->getText(id));
    std::string word;
    std::vector<T> tmp(values);
    while (multiValueStr>>word) {
        std::stringstream valueStr(word);
        T value;
        if (!(valueStr >> value)) {
            return false;
        }
        tmp.push_back(value);
    }
    values.swap(tmp);
    return true;
}

template <typename T, int N>
void XMLreaderProxy::read(Array<T,N>& values) const {
    if (!reader) return;
    std::stringstream multiValueStr(reader->getText(id));
    std::string word;
    values.resetToZero();
    int i=0;
    while (multiValueStr>>word && i<N) {
        std::stringstream valueStr(word);
        T value;
        if (!(valueStr >> value)) {
            EPC_ASSERT(false && (std::string("Cannot read value from XML element ") + reader->getName()).c_str());
        }
        values[i] = value;
        ++i;
    }
}

template <typename T, int N>
bool XMLreaderProxy::readNoThrow(Array<T,N>& values) const {
    if (!reader) return false;
    std::stringstream multiValueStr(reader->getText(id));
    std::string word;
    int i=0;
    while (multiValueStr>>word && i<N) {
        std::stringstream valueStr(word);
        T value;
        if (!(valueStr >> value)) {
            return false;
        }
        values[i] = value;
        ++i;
    }
    return true;
}

template<typename T> void XMLwriter::set(T const& value)
{
    std::stringstream valuestr;
    valuestr << value;
    valuestr >> data_map[currentId].text;
}

template<typename T> void XMLwriter::set(std::vector<T> const& values)
{
    std::stringstream valuestr;
    for (unsigned i=0; i<values.size(); ++i) {
        if (i != 0) {
            valuestr << " ";
        }
        valuestr << values[i];
    }
    data_map[currentId].text = valuestr.str();
}

template<typename T, int N> void XMLwriter::set(Array<T,N> const& values)
{
    std::stringstream valuestr;
    for (unsigned i=0; i<N; ++i) {
        if (i != 0) {
            valuestr << " ";
        }
        valuestr << values[i];
    }
    data_map[currentId].text = valuestr.str();
}

template<> inline void XMLwriter::set<bool>(bool const& value)
{
    if (value) {
        data_map[currentId].text = "True";
    }
    else {
        data_map[currentId].text = "False";
    }
}

template<typename ostrT>
void XMLwriter::toOutputStream(ostrT& ostr, int indent) const {
    if (data_map.empty()) return;

    if (isDocument) {
        ostr << "<?xml version=\"1.0\" ?>\n";
        std::vector<XMLwriter*> const& children = data_map.begin()->second.children;
        for (unsigned iNode=0; iNode<children.size(); ++iNode) {
            children[iNode]->toOutputStream(ostr);
        }
    }
    else {
        std::map<int,Data>::const_iterator it = data_map.begin();
        for (; it != data_map.end(); ++it) {
            std::vector<XMLwriter*> const& children = it->second.children;
            std::string const& text = it->second.text;
            std::string indentStr(indent, ' ');
            ostr << indentStr << "<" << name;
            if (data_map.size()>1 || it->first!=0) {
                ostr << " id=\"" << it->first << "\"";
            }
            if (children.empty()) {
                ostr << ">";
            }
            else {
                ostr << ">\n";
            }
            if (!text.empty()) {
                if (children.empty()) {
                    ostr << " " << text << " ";
                }
                else {
                    ostr << indentStr << "    " << text << "\n";
                }
            }
            for (unsigned iNode=0; iNode<children.size(); ++iNode) {
                children[iNode]->toOutputStream(ostr, indent+4);
            }
            if (!children.empty()) {
                ostr << indentStr;
            }
            ostr << "</" << name << ">\n";
        }
    }
}

}  // namespace plb

#endif  // XML_IO_HH
