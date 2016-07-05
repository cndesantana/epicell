#ifndef STATISTICS_HH
#define STATISTICS_HH

#include "statistics.h"
#include <cmath>

namespace epc {

namespace statistics {

Convergence::Convergence(double threshold_, long int maxSize_) 
	: threshold(threshold_), maxSize(maxSize_)
{ }

bool Convergence::hasConverged() {
	double avg = computeAverage();
	double std = computeStd();

	if ((long int)vec.size() < maxSize)
		return false;

	if (fabs(avg) > 1.0e-10){
		///std::cout << "vecSize = " << vec.size() << ", avg = " << avg << ", std = " << std << ", fabs(std/ avg) = " << fabs(std/avg) << " VS. threshold = " << threshold << std::endl;
		return fabs(std/ avg) < threshold;
	}
	else {
		//std::cout << "vecSize = " << vec.size() << ", avg = " << avg << ", std = " << std << ", fabs(std) = " << fabs(std) << " VS. threshold = " << threshold << std::endl;
		return fabs(std) < threshold;	
	}			
}

double Convergence::getNormalizedStd() {
	double avg = computeAverage();
	double std = computeStd();

	if (fabs(avg) > 1.0e-10){
		//std::cout << "vecSize = " << vec.size() << ", avg = " << avg << ", std = " << std << ", fabs(std/ avg) = " << fabs(std/avg) << " VS. threshold = " << threshold << std::endl;
		return fabs(std/ avg);
	}
	else {
		//std::cout << "vecSize = " << vec.size() << ", avg = " << avg << ", std = " << std << ", fabs(std) = " << fabs(std) << " VS. threshold = " << threshold << std::endl;
		return fabs(std);	
	}				
}

void Convergence::takeValue(double val) {
	if ((long int)vec.size() < maxSize) {
		vec.push_back(val);
	} else {
		vec.erase(vec.begin());
		vec.push_back(val);
	}
}

double Convergence::computeAverage() {
	double avg = double(0);

	for (auto v : vec) {
		avg += v;
	}
	return avg / (double)vec.size();
}

double Convergence::computeStd() {
	double avg = computeAverage();
	double std = (double)0;

	for (auto v : vec) {
		std += (v-avg)*(v-avg);
	}
	std /= (double)vec.size();

	return sqrt(std);
}


} // namespace statistics

} // namespace epc

#endif
