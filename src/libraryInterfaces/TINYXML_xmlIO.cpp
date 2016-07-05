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
 * Input/Output in XML format -- non-generic code.
 */

#include "core/epcDebug.h"
#include "mathematics/array.h"
#include "libraryInterfaces/TINYXML_xmlIO.h"
#include "libraryInterfaces/TINYXML_xmlIO.hh"
#include <fstream>

namespace epc {

XMLreader XMLreader::notFound;

XMLreader::XMLreader( std::vector<TiXmlNode*> pParentVect )
{
    mainProcessorIni(pParentVect);
}

XMLreader::XMLreader( std::string fName )
{
    TiXmlDocument* doc = 0;
    int loadOK = false;
    doc = new TiXmlDocument(fName.c_str());
    loadOK = doc->LoadFile();

    EPC_MSG_ASSERT(std::string("Problem processing input XML file ")+fName+". ", loadOK);

    mainProcessorIni(doc);
    delete doc;
}

void XMLreader::mainProcessorIni( TiXmlNode* pParent ) {
    std::vector<TiXmlNode*> pParentVect;
    pParentVect.push_back(pParent);
    mainProcessorIni(pParentVect);
}

void XMLreader::mainProcessorIni( std::vector<TiXmlNode*> pParentVect )
{
    std::map<int, TiXmlNode*> parents;
    for (unsigned iParent=0; iParent<pParentVect.size(); ++iParent) {
        EPC_ASSERT( pParentVect[iParent]->Type()==TiXmlNode::DOCUMENT ||
                          pParentVect[iParent]->Type()==TiXmlNode::ELEMENT );

        TiXmlElement* pParentElement = pParentVect[iParent]->ToElement();
        int id=0;
        if (pParentElement) {
            const char* attribute = pParentElement->Attribute("id");
            if (attribute) {
                std::stringstream attributestr(attribute);
                attributestr >> id;
            }
        }
        parents[id] = pParentVect[iParent];
    }

    std::map<int, TiXmlNode*>::iterator it = parents.begin();
    name = it->second->ValueStr();

    for (; it != parents.end(); ++it) {
        int id = it->first;

        TiXmlNode* pParent = it->second;
        Data& data = data_map[id];
        data.text="";

        typedef std::map<std::string, std::vector<TiXmlNode*> > ChildMap;
        ChildMap childMap;
        TiXmlNode* pChild;
        for ( pChild = pParent->FirstChild(); pChild != 0; pChild = pChild->NextSibling()) 
        {
            int type = pChild->Type();
            if ( type==TiXmlNode::ELEMENT ) {
                std::string name(pChild->Value());
                childMap[name].push_back(pChild);
            }
            else if ( type==TiXmlNode::TEXT ) {
                data.text = pChild->ToText()->ValueStr();
            }
        }
        for (ChildMap::iterator it = childMap.begin(); it != childMap.end(); ++it) {
            std::vector<TiXmlNode*> pChildVect = it->second;
            data.children.push_back( new XMLreader( pChildVect ) );
        }
    }
}

XMLreader::XMLreader() {
    name = "XML node not found";
}

XMLreader::~XMLreader() {
    std::map<int,Data>::iterator it = data_map.begin();
    for (; it != data_map.end(); ++it) {
        std::vector<XMLreader*>& children = it->second.children;
        for (unsigned iNode=0; iNode<children.size(); ++iNode) {
            delete children[iNode];
        }
    }
}

void XMLreader::print(int indent) const {
    std::string indentStr(indent, ' ');
    std::cout << indentStr << "[" << name << "]" << std::endl;
    std::string text = getFirstText();
    if (!text.empty()) {
        std::cout << indentStr << "  " << text << std::endl;
    }
    std::vector<XMLreader*> const& children = data_map.begin()->second.children;
    for (unsigned iNode=0; iNode<children.size(); ++iNode) {
        children[iNode]->print(indent+2);
    }
}

XMLreaderProxy XMLreader::operator[] (std::string name) const
{
    Data const& data = data_map.begin()->second;
    for (unsigned iNode=0; iNode<data.children.size(); ++iNode) {
        if (data.children[iNode]->name == name) {
            return XMLreaderProxy(data.children[iNode]);
        }
    }
    EPC_MSG_ASSERT(std::string("Element ")+name+std::string(" not found in XML file."), false);
    return XMLreaderProxy(0);
}

XMLreaderProxy XMLreader::getElement(std::string name, int id) const {
    std::map<int,Data>::const_iterator it = data_map.find(id);
    if (it==data_map.end()) {
        std::stringstream idStr;
        idStr << id;
        EPC_MSG_ASSERT(std::string("Element with id ") + idStr.str() + std::string(" does not exist"), false);
    }
    std::vector<XMLreader*> const& children = it->second.children;
    for (unsigned iNode=0; iNode<children.size(); ++iNode) {
        if (children[iNode]->name == name) {
            return XMLreaderProxy(children[iNode]);
        }
    }
    EPC_MSG_ASSERT(std::string("Element ")+name+std::string(" not found in XML file."), false);
    return XMLreaderProxy(0);
}

std::string XMLreader::getName() const {
    return name;
}

std::string XMLreader::getText() const {
    return data_map.begin()->second.text;
}

std::string XMLreader::getText(int id) const {
    std::map<int,Data>::const_iterator it = data_map.find(id);
    if (it != data_map.end()) {
        return it->second.text;
    }
    else {
        return "";
    }
}

int XMLreader::getFirstId() const {
    return data_map.begin()->first;
}

std::string XMLreader::getFirstText() const {
    return data_map.begin()->second.text;
}

bool XMLreader::idExists(int id) const {
    std::map<int,Data>::const_iterator it = data_map.find(id);
    if (it != data_map.end()) {
        return true;
    }
    else {
        return false;
    }
}

bool XMLreader::getNextId(int& id) const {
    std::map<int,Data>::const_iterator it = data_map.find(id);
    if (it != data_map.end()) {
        ++it;
        if (it != data_map.end()) {
            id = it->first;
            return true;
        }
    }
    return false;
}

std::vector<XMLreader*> const& XMLreader::getChildren(int id) const
{
    std::map<int,Data>::const_iterator it = data_map.find(id);
    if (it==data_map.end()) {
        EPC_MSG_ASSERT(std::string("Cannot access id ")+std::to_string(id)+" in XML element " + name, false);
    }
    return it->second.children;
}

XMLreaderProxy::XMLreaderProxy(XMLreader const* reader_)
    : reader(reader_)
{
    if (reader) {
        id = reader->getFirstId();
    }
    else {
        id = 0;
    }
}

XMLreaderProxy::XMLreaderProxy(XMLreader const* reader_, int id_)
    : reader(reader_),
      id(id_)
{ }

XMLreaderProxy XMLreaderProxy::operator[] (std::string name) const
{
    if (!reader) {
        EPC_MSG_ASSERT(std::string("Cannot read value from XML element ") + name, false);
    }
    return reader->getElement(name, id);
}

XMLreaderProxy XMLreaderProxy::operator[] (int newId) const {
    if (!reader) {
        EPC_ASSERT(false && "Cannot read value from XML element" );
    }
    if (!reader->idExists(newId)) {
        std::stringstream newIdStr;
        newIdStr << newId;
        EPC_MSG_ASSERT(std::string("Id ") + newIdStr.str() +
                std::string(" does not exist in XML element"), false);
    }
    return XMLreaderProxy(reader, newId);
}

bool XMLreaderProxy::isValid() const {
    return reader;
}

int XMLreaderProxy::getId() const {
    return id;
}

XMLreaderProxy XMLreaderProxy::iterId() const {
    if (!reader) {
        EPC_ASSERT(false && "Cannot read value from XML element" );
    }
    int newId = id;
    if (reader->getNextId(newId)) {
        return XMLreaderProxy(reader, newId);
    }
    else {
        return XMLreaderProxy(0);
    }
}

std::string XMLreaderProxy::getName() const {
    if (!reader) {
        EPC_ASSERT(false && "Cannot read value from XML element" );
    }
    return reader->getName();
}

std::vector<XMLreader*> const& XMLreaderProxy::getChildren() const {
    if (!reader) {
        EPC_ASSERT(false && "Cannot read value from XML element" );
    }
    return reader->getChildren(id);
}


XMLwriter::XMLwriter()
    : isDocument(true),
      currentId(0)
{ }

XMLwriter::~XMLwriter() {
    std::map<int,Data>::iterator it = data_map.begin();
    for (; it != data_map.end(); ++it) {
        std::vector<XMLwriter*>& children = it->second.children;
        for (unsigned iNode=0; iNode<children.size(); ++iNode) {
            delete children[iNode];
        }
    }
}

XMLwriter::XMLwriter(std::string name_)
    : isDocument(false),
      name(name_),
      currentId(0)
{ }

void XMLwriter::setString(std::string const& value)
{
    data_map[currentId].text = value;
}


XMLwriter& XMLwriter::operator[] (std::string name) {
    std::vector<XMLwriter*>& children = data_map[currentId].children;
    // If node already exists, simply return it.
    for (unsigned iNode=0; iNode<data_map[currentId].children.size(); ++iNode) {
        if (data_map[currentId].children[iNode]->name == name) {
            return *children[iNode];
        }
    }
    // Else, create and return it.
    children.push_back(new XMLwriter(name));
    return *children.back();
}

XMLwriter& XMLwriter::operator[] (int id) {
    currentId = id;
    return *this;
}

void XMLwriter::print(std::string fName) const {
    std::ofstream ofile(fName.c_str());
    EPC_MSG_ASSERT(std::string("Could not open file ") + fName
                                  + std::string(" for write access. "), !ofile.is_open());
    toOutputStream(ofile);
}


}  // namespace plb
