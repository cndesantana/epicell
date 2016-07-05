#ifndef EPITHELIUM_HH
#define EPITHELIUM_HH

#include "epithelium.h"
#include "cellCluster.h"

namespace epc {

// =============== Constuctors ========================= //
template<int d>
Epithelium<d>::Epithelium(CellCluster<d> &&cellCluster_) : cellCluster(cellCluster_)
{  }

// ========== Interface functions ================ //
template<int d>
CellCluster<d> const& Epithelium<d>::getCellCluster() const {
    return cellCluster;
}

template<int d>
std::vector<std::vector<std::shared_ptr<Vertex<d>>>> const& Epithelium<d>::getMovedVertices() const {
    return movedVertices;
}

template<int d>
CellCluster<d> & Epithelium<d>::getCellCluster() {
    return cellCluster;
}

template<int d>
std::vector<std::vector<std::shared_ptr<Vertex<d>>>> const& Epithelium<d>::getForcedVertices() const {
    return forcedVertices;
}

template<int d>
std::vector<std::shared_ptr<Vertex<d>>> const& Epithelium<d>::getMovedVertices(unsigned iA) const {
    return movedVertices[iA];
}

template<int d>
std::vector<std::shared_ptr<Vertex<d>>> const& Epithelium<d>::getForcedVertices(unsigned iA) const {
    return forcedVertices[iA];
}

template<int d>
template<class T, class U>
void Epithelium<d>::addSpecialVerticesAndForceFunction(T &&vertices, U &&foo) 
{
    forcedVertices.push_back(vertices);
    forceFunctions.push_back(foo);
}

template<int d>
template<class T, class U>
void Epithelium<d>::addSpecialVerticesAndMovementFunction(T &&vertices, U &&foo) 
{
    movedVertices.push_back(vertices);
    movementFunctions.push_back(foo);
}

template<int d>
template<class T>
void Epithelium<d>::addFixedOnYVertices(T &&vertices) 
{
    for (auto v: vertices){
        fixedOnYVertices.push_back(v);
    }
}

template<int d>
void Epithelium<d>::removeAllForcedVerticesAndForceFunction() 
{
    forcedVertices.clear();
    forceFunctions.clear();
}

template<int d>
void Epithelium<d>::removeAllMovedVerticesAndMovementFunction() 
{
    movedVertices.clear();
    movementFunctions.clear();
}

template<int d>
void Epithelium<d>::removeAllFixedOnYVertices() 
{
    fixedOnYVertices.clear();
}

// ========== dynamic functions ================ //
template<int d>
void Epithelium<d>::performTimeStep(double dt) {
    moveSpecialVertices(dt);
    cellCluster.resetForce();
    cellCluster.updateForce();
    addForceOnSpecialVertices();
    setNullForceOnMovedVertices();
    setNullForceYOnFixedYVertices();
    cellCluster.advance(dt);
}

template<int d>
void Epithelium<d>::updateCellMechanics(double gamma, double lambda, double sdGamma, double sdLambda){
    for(auto c: cellCluster.getCells()){
        c->setGamma(gamma,sdGamma);
    }
    for(auto e: cellCluster.getEdges()){
        e->setLambda(lambda,sdLambda);
    }
    return; 
}

template<int d>
void Epithelium<d>::updateForce() {
    cellCluster.resetForce();
    cellCluster.updateForce();
}
// ========== dynamic functions ================ //
template<int d>
void Epithelium<d>::addForceOnSpecialVertices() {
    EPC_ASSERT(forceFunctions.size() ==  forcedVertices.size());
    for (unsigned iA = 0; iA < forcedVertices.size(); ++iA) {
        for (auto v : forcedVertices[iA]) {
            v->addToForce(forceFunctions[iA](*v));
        }
    }
}

template<int d>
void Epithelium<d>::moveSpecialVertices(double dt) {
    EPC_ASSERT(movementFunctions.size() ==  movedVertices.size());
    for (unsigned iA = 0; iA < movedVertices.size(); ++iA) {
        for (auto v : movedVertices[iA]) {
            v->setPosition(v->getPosition()+movementFunctions[iA](*v)*dt);
        }

    }
}

template<int d>
void Epithelium<d>::setNullForceOnMovedVertices() {
    Array<double,d> zero(0.0);
    for (unsigned iA = 0; iA < movedVertices.size(); ++iA) {
        for (auto v : movedVertices[iA]) {
            v->setForce(zero);
        }

    }
}

template<int d>
void Epithelium<d>::setNullForceYOnFixedYVertices() {
    Array<double,d> nullForceY(0.0);
    for (auto v: fixedOnYVertices) {
        nullForceY[0] = v->getForce()[0];
        v->setForce(nullForceY);
    }
}

}  // end namespace epc

#endif
