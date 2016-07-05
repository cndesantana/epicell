#ifndef VERTEX_HH
#define VERTEX_HH

#include <memory>
#include "mathematics/array.h"
#include "vertex.h"

namespace epc {

// =============== Constuctors ========================= //
template<int d>
Vertex<d>::Vertex() : position(Array<double,d>(double())), position_t1(position), 
					  force(Array<double,d>(double())), dyn(0)
{ 
    createNewId();
}

template<int d>
Vertex<d>::Vertex(Array<double,d> position_) : 
    position(position_), position_t1(position), force(Array<double,d>(double())), dyn(0)
{ 
    createNewId();
}

template<int d>
Vertex<d>::Vertex(Array<double,d> position_, std::shared_ptr<Dynamics<d> > dyn_) : 
	position(position_), position_t1(position), force(Array<double,d>(double())), dyn(dyn_)
{ 
    createNewId();
}

template<int d>
Vertex<d>::Vertex(Array<double,d> position_, Array<double,d> force_, std::shared_ptr<Dynamics<d> > dyn_) : 
	position(position_), position_t1(position), force(force_), dyn(dyn_)
{ 
    createNewId();
}

template<int d>
void Vertex<d>::createNewId() {
    static unsigned long allIds = 0; 
    id = allIds++;
}

// ========== Interface functions ================ 
template<int d>
Dynamics<d> const& Vertex<d>::getDynamics() const {
	EPC_ASSERT(dyn);
	return *dyn;
}

template<int d>
Dynamics<d> & Vertex<d>::getDynamics() {
	EPC_ASSERT(dyn);
	return *dyn;
}

template<int d>
Array<double,d> const& Vertex<d>::getForce() const {
	return force;
}

template<int d>
Array<double,d> const& Vertex<d>::getPosition() const {
	return position;
}

template<int d>
unsigned long const& Vertex<d>::getId() const {
	return id;
}

// ========== update functions =================== //
template<int d>
void Vertex<d>::addToForce(const Array<double,d> &force_) {
	force += force_;
}

template<int d>
void Vertex<d>::setForce(const Array<double,d> &force_) {
	force = force_;
}

template<int d>
void Vertex<d>::setPosition(const Array<double,d> &position_) {
	position = position_;
}

template<int d>
void Vertex<d>::setDynamics(const std::shared_ptr<Dynamics<d> > dyn_) {
	dyn = dyn_;
}

// ========== Kinematic functions ================ //
template<int d>
void Vertex<d>::advance(double dt) {
	
	dyn->advance(position,position_t1,force,dt); // here the dynamics will divide the force by the mass to get the acceleration.
}

}  // end namespace epc

#endif
