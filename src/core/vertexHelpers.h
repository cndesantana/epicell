#ifndef VERTEX_HELPERS_H
#define VERTEX_HELPERS_H

#include <memory>
#include "mathematics/array.h"
#include "mathematics/constants.h"
#include <algorithm>    // std::sort
#include <math.h>

namespace epc {

template<unsigned d>
struct SortHelpers {
    SortHelpers(Array<double,d> centroid_) : centroid(centroid_)
    { }

    bool operator() (const std::shared_ptr<Vertex<d>> &v1,const std::shared_ptr<Vertex<d>> &v2) { 
        Array<double,d> p1 = v1->getPosition()-centroid;
        Array<double,d> p2 = v2->getPosition()-centroid;

        double theta1 = atan2(p1[1],p1[0]);
        double theta2 = atan2(p2[1],p2[0]);
        return theta1 < theta2;
    }
private:
    Array<double,d> centroid;
};

template<unsigned d>
struct vertexHelpers {
static void sortVerticesCounterClockwiseInPlace(std::vector<std::shared_ptr<Vertex<d>>> &vertices) {
    static_assert(d == 2, "sorting not implemented in 3d yet.");
    Array<double,d> centroid = computeCentroid(vertices);
    std::sort(vertices.begin(),vertices.end(), SortHelpers<d>(centroid));
}

static std::vector<std::shared_ptr<Vertex<d>>> sortVerticesCounterClockwise(const std::vector<std::shared_ptr<Vertex<d>>> &vertices) {
    static_assert(d == 2, "sorting not implemented in 3d yet.");

    std::vector<std::shared_ptr<Vertex<d>>> cpVertices = vertices;
    sortVerticesCounterClockwiseInPlace(cpVertices);
    return cpVertices;
}

static Array<double,d> computeCentroid(const std::vector<std::shared_ptr<Vertex<d>>> &vertices) {
    Array<double,d> centroid(double(0.0));
    for ( auto vertex : vertices ) {
        centroid += vertex->getPosition();
    }
    centroid /= (double)vertices.size();
    return centroid;
}

static std::vector<std::shared_ptr<Vertex<d>>> generateHexagon(const Array<double,d> &center, double radius) {
    static_assert(d == 2, "d = 3 not implemented yed");
    static std::vector<std::shared_ptr<Vertex<d>>> vertices;
    for (auto theta : mathConstants::HexagonAngles) {
        Array<double,d> pos;
        pos[0] = radius * cos(theta);
        pos[1] = radius * sin(theta);
        vertices.push_back(std::shared_ptr<Vertex<d>>(new Vertex<d>(pos)));
    }
    return vertices;
}

static bool isInVector(Vertex<d> const &v, const std::vector<std::shared_ptr<Vertex<d>>> &vertices) {
    for (auto vTmp : vertices) {
        if (v.getId() == vTmp->getId()) return true;
    }
    return false;
}

static void appendIfNotInVector(const std::shared_ptr<Vertex<d>> &v, std::vector<std::shared_ptr<Vertex<d>>> &vertices) {
    if (!isInVector(*v, vertices)) vertices.push_back(v);
}

/// Checks the equality of position of two vertices
static bool isEqual(Vertex<d> const &v1, Vertex<d> const &v2, 
    double threshold = (double)100*std::numeric_limits<double>::epsilon())
{
    return arrayOperations<double,d>::isEqual(v1.getPosition(),v2.getPosition(),threshold);
}

}; // end vertex helper

}  // end namespace epc

#endif
