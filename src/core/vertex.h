#ifndef VERTEX_H
#define VERTEX_H

#include <memory>
#include "mathematics/array.h"
#include "dynamics.h"

namespace epc {

template<int d>
class Vertex {
public:
	/// Default constructor (everything is set to zero)
	Vertex();
	Vertex(Array<double,d> position_, std::shared_ptr<Dynamics<d> > dyn_);
	Vertex(Array<double,d> position_);
	Vertex(Array<double,d> position_, Array<double,d> force_, std::shared_ptr<Dynamics<d> > dyn_);
private:
    static_assert(d <= 3 && d >= 2, "d must be 2 or 3!");
    void createNewId();
public: 
// ========== Interface functions ================ //

    /// Read-write access to vertex positions.
    /** \param iD index of the accessed distribution function */
    double& operator[](unsigned iD) {
        return position[iD];
    }
    /// Read-only access to vertex positions.
    /** \param iD index of the accessed distribution function */
    double const& operator[](unsigned iD) const {
    	return position[iD];
    }

    /// Get a reference to non-modifiable dynamics
    Dynamics<d> const& getDynamics() const;

    /// Get a reference to modifiable dynamics
    Dynamics<d>& getDynamics();

    /// Get a reference to non-modifiable force
    Array<double,d> const& getForce() const;

    /// Get a reference to non-modifiable position
    Array<double,d> const& getPosition() const;

    /// Get a reference to non-modifiable id
    unsigned long const& getId() const;

// ========== update functions =================== //
    /// adding a force to the force acting on the vertex
    void addToForce(const Array<double,d> &force_);
    /// sets a new force to the vertex
    void setForce(const Array<double,d> &force_);

    /// sets a new position to the vertex
    void setPosition(const Array<double,d> &position_);
    /// sets a new dynamics
    void setDynamics(const std::shared_ptr<Dynamics<d> > dyn_);

// ========== Kinematic functions ================ //
	/// executes a time-step for the vertex. 
	/// modifies the position of the vertex.
	void advance(double dt);

private:
	Array<double,d> position, position_t1, force; // position and force acting on each vertex
	std::shared_ptr<Dynamics<d> > dyn;	   // pointer on dynamics object which handles the movement of each vertex
    unsigned long id; // unique identifier for each Edge TODO: rethink that maybe
};

}  // end namespace epc

#endif
