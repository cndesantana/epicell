#ifndef DYNAMICS_HH
#define DYNAMICS_HH

#include "dynamics.h"

namespace epc {

template <int d>
void NoDynamics<d>::advance(Array<double,d> &position, Array<double,d> &position_t1, const Array<double,d> &force, double dt)
{ 
	// Does nothing.
}

template <int d>
Verlet<d>::Verlet(double damping_, double mass_) : Dynamics<d>(), damping(damping_), mass(mass_)
{ 
	EPC_ASSERT(damping >= double(0.0) && "damping must be higher or equal to 0.")
	EPC_ASSERT(mass > double(0.0) && "mass must be higher than 0.")
}

template <int d>
void Verlet<d>::advance(Array<double,d> &position, 
		Array<double,d> &position_t1, const Array<double,d> &force, double dt)
{
	Array<double,d> acceleration = force/mass;
	Array<double,d> old_pos = position;
	position = (double(2)-damping*dt)*position - (double(1)-damping*dt)*position_t1 + dt*dt*acceleration;
	position_t1 = old_pos;
}

}  // end namespace epc

#endif