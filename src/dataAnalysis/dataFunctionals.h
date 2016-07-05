#ifndef DATA_FUNCTIONALS_H
#define DATA_FUNCTIONALS_H

#include "mathematics/constants.h"
#include "core/cellCluster.h"
#include "environment/elasticEnvironment.h"
#include <cmath>

namespace epc {

namespace dataFunctionals2D {
    
// ============== cellCluster functionals ================== //
std::function<std::vector<double> (const CellCluster<2> &cellCluster)> computeArea = [] (const CellCluster<2> &cellCluster){
    std::vector<double> vals;
    for (auto c: cellCluster.getCells()) {
        const auto verts = vertexHelpers<2>::sortVerticesCounterClockwise(cellHelpers<2>::getVertices(*c));
        for (unsigned iA = 0; iA < verts.size(); ++iA) {
            vals.push_back(c->computeArea());
        }
    }
    return vals;
};


std::function<std::vector<double> (const CellCluster<2> &cellCluster)> computeDegree = [] (const CellCluster<2> &cellCluster){
    std::vector<double> vals;
    for (auto c: cellCluster.getCells()) {
        vals.push_back(c->computeDegree());
    }
    return vals;
};


std::function<std::vector<double> (const CellCluster<2> &cellCluster)> isOnMitosis = [] (const CellCluster<2> &cellCluster){
    std::vector<double> vals;
    for (auto c: cellCluster.getCells()) {
        vals.push_back(c->isOnMitosis());
    }
    return vals;
};

std::function<std::vector<double> (const CellCluster<2> &cellCluster)> getType = [] (const CellCluster<2> &cellCluster){
    std::vector<double> vals;
    for (auto c: cellCluster.getCells()) {
        vals.push_back(c->getType().getId());
    }
    return vals;
};

std::function<std::vector<double> (const CellCluster<2> &cellCluster)> getProlifSignal = [] (const CellCluster<2> &cellCluster){
    std::vector<double> vals;
    for (auto c: cellCluster.getCells()) {
        vals.push_back(c->getSignal());
    }
    return vals;
};

std::function<std::vector<double> (const CellCluster<2> &cellCluster)> getSqrProlifSignal = [] (const CellCluster<2> &cellCluster){
    std::vector<double> vals;
    for (auto c: cellCluster.getCells()) {
        vals.push_back(pow(c->getSignal(),2));
    }
    return vals;
};
std::function<std::vector<double> (const CellCluster<2> &cellCluster)> getCubicProlifSignal = [] (const CellCluster<2> &cellCluster){
    std::vector<double> vals;
    for (auto c: cellCluster.getCells()) {
        vals.push_back(pow(c->getSignal(),3));
    }
    return vals;
};

std::function<std::vector<Array<double,2>> (const CellCluster<2> &cellCluster)> computeForce = [] (const CellCluster<2> &cellCluster) {
    std::vector<Array<double,2>> vals;
    for (auto c: cellCluster.getCells()) {
        const auto verts = vertexHelpers<2>::sortVerticesCounterClockwise(cellHelpers<2>::getVertices(*c));
        for (auto v: verts) {
            vals.push_back(v->getForce());
        }
    }
    return vals;
};

// ============== force functionals ================== //

auto constForce = [](const Vertex<2> &v, const Array<double,2> &f) -> Array<double,2> { return f;};

// ============== force functionals ================== //

auto radialForce = [](const Vertex<2> &v, const Array<double,2> &center, double amplitude) -> Array<double,2> 
{ 
    Array<double,2> r = v.getPosition()-center;
    double distance = r.norm();
    r /= distance;

    return amplitude * r;
};

// ============== force functionals ================== //

auto elasticShellForce = [](const Vertex<2> &v, const ElasticShellEnvironment<2> &elasticShell) -> Array<double,2> 
{ 
    Array<double,2> r = elasticShell.computeForceOnVertex(v);

    return r;
};

} // end dataFunctionals2D


}  // end namespace epc

#endif
