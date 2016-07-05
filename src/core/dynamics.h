#ifndef DYNAMICS_H
#define DYNAMICS_H

#include "mathematics/array.h"

namespace epc {

template <int d>
class Dynamics {
public:
	/// Destructor: virtual to enable inheritance
    virtual ~Dynamics() { }
    virtual void advance(Array<double,d> &position, Array<double,d> &position_t1, const Array<double,d> &force, double dt) = 0;
};

template <int d>
class NoDynamics final : public Dynamics<d> {
	virtual void advance(Array<double,d> &position, Array<double,d> &position_t1, const Array<double,d> &force, double dt) override;
};

template <int d>
class Verlet final : public Dynamics<d> {
public:
	Verlet(double damping_ = double(0.0), double mass_ = double(1.0));

	virtual void advance(Array<double,d> &position, 
		Array<double,d> &position_t1, const Array<double,d> &force, double dt) override;

private:
	double damping;
	double mass;
};


}  // end namespace epc

#endif
