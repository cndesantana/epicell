#ifndef EDGE_HH
#define EDGE_HH

#include "vertex.h"
#include "edge.h"

namespace epc {

// =============== Constuctors ========================= //
template<int d>
Edge<d>::Edge() : vertices(std::tuple<std::shared_ptr<Vertex<d> > , std::shared_ptr<Vertex<d> > >(0,0)), 
                  lambda(double(), sdLambda(double()))
{
    createNewId();
}

template<int d>
Edge<d>::Edge(std::shared_ptr<Vertex<d> > &vertexOne_, std::shared_ptr<Vertex<d> > &vertexTwo_, double lambda_,double sdLambda_) : 
        vertices(std::make_tuple(vertexOne_,vertexTwo_)), 
        lambda(lambda_),
        sdLambda(sdLambda_)
{
    EPC_ASSERT(getVertexOne()->getId() !=  getVertexTwo()->getId() && "An edge must be made of two different vertices");
    EPC_ASSERT(getVertexOne()->getPosition() !=  getVertexTwo()->getPosition() && "An edge must be made of vertices with different positions");
    createNewId();
    setLambda(lambda_, sdLambda_);
}

template<int d>
void Edge<d>::createNewId() {
    static unsigned long allIds = 0; 
    id = allIds++;
}

// ========== Interface functions ================ //
template<int d>
unsigned long const& Edge<d>::getId() const {
    return id;
}

template<int d>
std::shared_ptr<Vertex<d>> const& Edge<d>::getVertexOne() const {
    EPC_ASSERT(std::get<0>(vertices));
    return (std::get<0>(vertices));
}

template<int d>
std::shared_ptr<Vertex<d>> & Edge<d>::getVertexOne() {
    EPC_ASSERT(std::get<0>(vertices));
    return (std::get<0>(vertices));
}

template<int d>
std::shared_ptr<Vertex<d>> const& Edge<d>::getVertexTwo() const {
    EPC_ASSERT(std::get<1>(vertices));
    return (std::get<1>(vertices));
}

template<int d>
std::shared_ptr<Vertex<d>> & Edge<d>::getVertexTwo() {
    EPC_ASSERT(std::get<1>(vertices));
    return (std::get<1>(vertices));
}

template<int d>
double const& Edge<d>::getLambda() const {
    return lambda;
}

template<int d>
double const& Edge<d>::getSdLambda() const {
    return sdLambda;
}

template<int d>
void Edge<d>::setLambda(double lambda_,double sdLambda_) {
//    lambda = lambda_;
    sdLambda = sdLambda_;
    lambda = NRNG(lambda_,sdLambda_).ngenerate();
}

template<int d>
void Edge<d>::setId(unsigned long id_) {
    id = id_;
}

// ========== geometric functions ================ //
template<int d>
double Edge<d>::computeLength() {
    if (computeVectorFromOneToTwo().norm() != computeVectorFromTwoToOne().norm()){
        std::cout << "ERROR: Edge " << getId() << std::endl;
        std::cout << "VertexOne = " << getVertexOne()->getPosition() << " vs. VertexTwo = " << getVertexTwo()->getPosition() << std::endl;
        std::cout << "computeVectorFromOneToTwo() = " << computeVectorFromOneToTwo() << " vs. computeVectorFromTwoToOne() = " << computeVectorFromTwoToOne() << std::endl;
        std::cout << "computeVectorFromOneToTwo().norm() = " << computeVectorFromOneToTwo().norm() << " vs. computeVectorFromTwoToOne().norm() = " << computeVectorFromTwoToOne().norm() << std::endl;
    }
    EPC_ASSERT(computeVectorFromOneToTwo().norm() == computeVectorFromTwoToOne().norm());
    return computeVectorFromOneToTwo().norm();
}

template<int d>
Array<double,d> Edge<d>::computeVectorFromOneToTwo() {
    return getVertexTwo()->getPosition()-getVertexOne()->getPosition();
}

template<int d>
Array<double,d> Edge<d>::computeVectorFromTwoToOne() {
    return getVertexOne()->getPosition()-getVertexTwo()->getPosition();
}

// ========== dynamic functions ================ //
template<int d>
void Edge<d>::computeAndUpdateForce() {

    Array<double,d> force(computeVectorFromOneToTwo());
    force /= force.norm();
    force *= lambda;
    getVertexOne()->addToForce(force);
    getVertexTwo()->addToForce(-force);    
}

}  // end namespace epc

#endif
