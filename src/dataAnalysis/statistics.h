#ifndef STATISTICS_H
#define STATISTICS_H

#include <vector>
#include <cmath>

namespace epc {

namespace statistics {

class Convergence {
public:
	Convergence(double threshold_, long int maxSize_);

public:
	bool hasConverged();
	void takeValue(double val);
	double getNormalizedStd();
private:
	double computeAverage();
	double computeStd();
private:
	double threshold;
	long int maxSize;
	std::vector<double> vec;
};


} // namespace statistics

} // namespace epc

#endif
