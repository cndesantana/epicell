#ifndef SIGNALING_H
#define SIGNALING_H


#include "core/cellCluster.h"

namespace epc {

template <int d>
class Signaling {
public:
    Signaling(double dissipCoef_= double(0.0), double propagCoef_=double(0.0));
    void propagateSignal(CellCluster<d> &cellCluster);
private:

	double const& getDissipCoef() const;
	double const& getPropagCoef() const;

	void setDissipCoef(double dissipCoef_);
	void setPropagCoef(double propagCoef_);

	void setNullReceivedSignal(CellCluster<d> &cellCluster);
	void computePropagatedSignal(CellCluster<d> &cellCluster);
	void updateSignal(CellCluster<d> &cellCluster);

private:
	double dissipCoef;
	double propagCoef;
};


}  // end namespace epc

#endif
