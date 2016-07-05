#ifndef IO_H
#define IO_H

#include <iostream>
#include "mathematics/array.h"
#include "core/cell.h"
#include "core/epcDebug.h"

namespace epc {

template<typename T, unsigned size>
std::ostream& operator<< (std::ostream &out, const Array<T,size> &ary) {
    for (unsigned iD = 0; iD < size-1; ++iD) {
        out << ary[iD] << " ";
    }
    out << ary[size-1];
    return out;
}

template<int d>
std::ostream& operator<< (std::ostream &out, const Cell<d> &cell) {
    const auto verts = vertexHelpers<d>::sortVerticesCounterClockwise(cellHelpers<d>::getVertices(cell));

    for (auto v : verts) {
        out << v->getPosition() << std::endl;
    }
    return out;
}

} // end namespace epc

#endif

