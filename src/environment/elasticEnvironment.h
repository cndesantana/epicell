#ifndef ELASTICENVIRONMENT_H
#define ELASTICENVIRONMENT_H

#include <core/vertex.h>
#include <mathematics/array.h>

namespace epc {

template<int d>
class ElasticEnvironment {

public:
	ElasticEnvironment(double elasticCoef_);

	Array<double, d> computeForceOnVertex(const Vertex<d> &v) const;

private:
	virtual bool isPositionOutside(Array<double, d> position) const =0;
	virtual double computeDistanceFromBorders(Array<double, d> position) const =0 ;
	virtual Array<double, d> computeForceDirection(Array<double, d> position) const =0 ;

private:
	double elasticCoef;
};

template<int d>
class ElasticShellEnvironment : public ElasticEnvironment<d>{

public:
	ElasticShellEnvironment(double elasticCoef_, Array<double,2> center_, double radius_);
private:
	bool isPositionOutside(Array<double, d> position) const;
	double computeDistanceFromBorders(Array<double, d> position) const;
	Array<double, d> computeForceDirection(Array<double, d> position) const;

private:
	Array<double, d> center;
	double radius;
};
}  // end namespace epc

#endif
