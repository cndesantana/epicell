#ifndef ELASTICENVIRONMENT_HH
#define CELLTYPE_HH

#include "elasticEnvironment.h"

namespace epc {

template<int d>
ElasticEnvironment<d>::ElasticEnvironment(double elasticCoef_) : elasticCoef(elasticCoef_){}

template<int d>
Array<double, d> ElasticEnvironment<d>::computeForceOnVertex(const Vertex<d> &v) const{
	// std::cout << "computeForceOnVertex?" << std::endl;
	if (isPositionOutside(v.getPosition())){
		return elasticCoef*computeDistanceFromBorders(v.getPosition())*computeForceDirection(v.getPosition());
	}
	return Array<double, d>(double(0.0));
	
}


template<int d>
ElasticShellEnvironment<d>::ElasticShellEnvironment(double elasticCoef_, Array<double,2> center_, double radius_) : ElasticEnvironment<d>(elasticCoef_), center(center_), radius(radius_){}



template<int d>
bool ElasticShellEnvironment<d>::isPositionOutside(Array<double, d> position) const{
	if ((position - center).norm() > radius){
		// std::cout << "isPositionOutside?" << std::endl;
		// std::cout << "\t YES!" << std::endl;
		return true;
	}
	// std::cout << "\t NO!" << std::endl;
	return false;
}

template<int d>
double ElasticShellEnvironment<d>::computeDistanceFromBorders(Array<double, d> position) const{
	// std::cout << "computeDistanceFromBorders?" << std::endl;
	return (position - center).norm();
}

template<int d>
Array<double, d> ElasticShellEnvironment<d>::computeForceDirection(Array<double, d> position) const{
	// std::cout << "computeForceDirection?" << std::endl;
	return (center - position)/((position - center).norm());
}

}  // end namespace epc

#endif
