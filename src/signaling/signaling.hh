#ifndef SIGNALING_HH
#define SIGNALING_HH


#include "signaling.h"

namespace epc {

// =============== Constuctors ========================= //
template<int d>
Signaling<d>::Signaling(double dissipCoef_, double propagCoef_) : dissipCoef(dissipCoef_), propagCoef(propagCoef_)
{}

// ========== Interface functions ================ 
template<int d>
double const& Signaling<d>::getDissipCoef() const {
	return dissipCoef;
}

template<int d>
double const& Signaling<d>::getPropagCoef() const {
	return propagCoef;
}

// ========== update functions =================== //

template<int d>
void Signaling<d>::setDissipCoef(double dissipCoef_) {
	dissipCoef = dissipCoef_;
}

template<int d>
void Signaling<d>::setPropagCoef(double propagCoef_) {
	propagCoef = propagCoef_;
}

template<int d>
void Signaling<d>::setNullReceivedSignal(CellCluster<d> &cellCluster) {
	for (auto c:cellCluster.getCells()){	
		c->setReceivedSignal(double(0.0));
	}
}

template<int d>
void Signaling<d>::propagateSignal(CellCluster<d> &cellCluster) {
	setNullReceivedSignal(cellCluster);
	computePropagatedSignal(cellCluster);
	updateSignal(cellCluster);
}

template<int d>
void Signaling<d>::computePropagatedSignal(CellCluster<d> &cellCluster) {
	auto cells = cellCluster.getCells();	
	auto connect = cellCluster.getConnectivity();
	for (unsigned iC = 0; iC < connect.size(); iC++){
		double cPerimeter = cells[iC]->computePerimeter();
		for (unsigned iA = 0; iA < connect[iC].size(); iA++){
			auto neighbor = std::get<0>(connect[iC][iA]);
			auto edge = std::get<1>(connect[iC][iA]);
			neighbor->addToReceivedSignal(propagCoef*cells[iC]->getSignal()*edge->computeLength()/cPerimeter);	
		}
	}
	
}

template<int d>
void Signaling<d>::updateSignal(CellCluster<d> &cellCluster) {
	for (auto c:cellCluster.getCells()){	
		c->setOldSignal(c->getSignal());
		c->setSignal((1-dissipCoef)*((1-propagCoef)*c->getSignal()+c->getReceivedSignal()));
	}
}


}  // end namespace epc

#endif
